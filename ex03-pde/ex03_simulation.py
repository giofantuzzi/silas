import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt
from scipy.io import savemat
import h5py

from reacDiff1D.helpers import (
    make_Tn_fields, 
    cheb_coeffs, 
    RunningMeans, 
    weighted_pod, 
    compute_quadrature_weights
)


transient = False  # True for running up to T_final
DNS = True        # True for solving PDE from T_final long enough to collect N snapshots
POD = True         # True to compute POD modes and save in .mat format 

Nx = 128
D = 4.0 * 0.0322307
eps = 0.01
alpha = 0.01
dealias = 2
timestepper = d3.RK443
timestep = 1e-3
dtype = np.float64

T_final = 100              
save_every = 5 #25
log_every = 5000
cheb_orders = list(range(3))     # n = 0..4

n_pod_modes = [3, 4, 5]

N_train = 25000 #250000
N_test = 20000 
N_total = N_train + N_test

data_dir = "./paper-code/reacDiff-pde/data"
plot_dir = "./paper-code/reacDiff-pde/results"
final_state_file = f"{data_dir}/final_state_Nx={Nx}_T={T_final:.2f}.npz"
means_plot_file  = f"{data_dir}/mean_an_T={T_final:.2f}.png"

out_npz = f"{data_dir}/reacDiff_full_N={N_total}.npz"

# If restarting: load from a saved final state file instead of using constant ICs
restart_from_file = None  # e.g. "./reacDiff1D/data/final_state_Nx={Nx}_T={T_final:.2f}.npz"


plot_singular_values = False


# -----------------------
# Main
# -----------------------
if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    xcoord = d3.Coordinate('x')
    dist = d3.Distributor(xcoord, dtype=dtype)
    xbasis = d3.Chebyshev(xcoord, size=Nx, bounds=(-1, 1), dealias=dealias)

    u = dist.Field(name='u', bases=xbasis)
    v = dist.Field(name='v', bases=xbasis)
    tau_u1 = dist.Field(name='tau_u1')
    tau_u2 = dist.Field(name='tau_u2')
    tau_v1 = dist.Field(name='tau_v1')
    tau_v2 = dist.Field(name='tau_v2')

    dx = lambda A: d3.Differentiate(A, xcoord)

    tau_basis = xbasis.derivative_basis(2)
    p1 = dist.Field(bases=tau_basis)
    p2 = dist.Field(bases=tau_basis)
    p1['c'][-1] = 1
    p2['c'][-2] = 2
    

    problem = d3.IVP([u, v, tau_u1, tau_u2, tau_v1, tau_v2], namespace=locals())
    problem.add_equation("dt(u) - D*dx(dx(u)) + tau_u1*p1 + tau_u2*p2 = (v - u*u - u*u*u)/eps")
    problem.add_equation("dt(v) - D*dx(dx(v)) + tau_v1*p1 + tau_v2*p2 = -u + alpha")

    problem.add_equation("u(x='left') = -2")
    problem.add_equation("u(x='right') = -2")
    problem.add_equation("v(x='left') = -4")
    problem.add_equation("v(x='right') = -4")

    solver = problem.build_solver(timestepper)
    x = dist.local_grid(xbasis)

    if transient:
        # Initial condition
        if restart_from_file is None:
            u['g'] = -2
            v['g'] = -4
            logger.info("Using constant initial condition.")
        else:
            data = np.load(restart_from_file)
            u['g'] = data["u_final"].copy()
            v['g'] = data["v_final"].copy()
            solver.sim_time = float(data["t_final"])
            logger.info(f"Restarting from {restart_from_file} at t={solver.sim_time:g}.")

        # Monitoring setup
        Tn_fields = make_Tn_fields(dist, xbasis, x, cheb_orders)
        running = RunningMeans(m=len(cheb_orders))

        # Run to target time
        t_start = solver.sim_time
        t_target = t_start + T_final
        logger.info(f"Running from t={t_start:g} to t={t_target:g} with dt={timestep:g}.")

        # Record at t_start
        a0 = cheb_coeffs(u, Tn_fields)
        running.update(solver.sim_time, a0)

        while solver.proceed and solver.sim_time < t_target - 0.5*timestep:
            solver.step(timestep)

            if (solver.iteration % log_every) == 0:
                logger.info(f"Iter={solver.iteration}, t={solver.sim_time:.6e}, ||u||={np.linalg.norm(u['g']):.3e}")

            if (solver.iteration % save_every) == 0:
                a_vec = cheb_coeffs(u, Tn_fields)
                running.update(solver.sim_time, a_vec)

        # Extract arrays
        t_arr, a_arr, mean_arr = running.arrays()
        t_arr = t_arr[1:]
        mean_arr = mean_arr[1:, :]
        # Plot mean_a_n(t)
        plt.figure()
        for j, n in enumerate(cheb_orders):
            plt.plot(t_arr, mean_arr[:, j], label=f"n={n}")
        plt.xlabel("t")
        plt.ylabel("mean_a_n(t)")
        plt.title("Cumulative time averages of Chebyshev coefficients")
        plt.legend()
        plt.tight_layout()
        # plt.savefig(means_plot_file, dpi=200)
        # logger.info(f"Saved plot to: {means_plot_file}")
        # plt.show()

        # Save final state for restart

        u.change_scales(1) 
        v.change_scales(1) 

        np.savez(
            final_state_file,
            x=x.copy(),
            t_final=np.array([solver.sim_time]).item(),
            u_final=u['g'].copy(),
            v_final=v['g'].copy(),
            params=np.array([Nx, D, eps, alpha, timestep], dtype=float),
            cheb_orders=np.array(cheb_orders, dtype=int),
        )
        logger.info(f"Saved final state to: {final_state_file}")

    if DNS:
        # restart from file with removed transient
        data = np.load(final_state_file)

        u['g'] = data["u_final"].copy()
        v['g'] = data["v_final"].copy()
        solver.sim_time = float(data["t_final"])
        # logger.info(f"Restarted from {restart_from_file} at t={solver.sim_time:g}")
        logger.info(f"Restarted from t={solver.sim_time:g}")

        # ---- run long enough to collect N_total snapshots ----
        T_total = (N_total - 1) * save_every * timestep
        t_stop = solver.sim_time + T_total
        solver.stop_sim_time = t_stop

        u_list = [u['g'].copy()]
        v_list = [v['g'].copy()]
        t_list = [solver.sim_time]
        E_list = [d3.Integrate(0.5*(u*u + v*v), 'x').evaluate()['g'].item()]        
        
        u_rhs_list = []
        v_rhs_list = []

        # compute u_rhs and v_rhs at current time (t0) and store them
        dxu = d3.Differentiate(u, xcoord).evaluate()
        dxdxu = d3.Differentiate(dxu, xcoord).evaluate()
        u_rhs = ((v - u*u - u*u*u)/eps + D*dxdxu).evaluate()  
        u_rhs.change_scales(1)
        u_rhs_list.append(u_rhs['g'].copy())

        dxv = d3.Differentiate(v, xcoord).evaluate()
        dxdxv = d3.Differentiate(dxv, xcoord).evaluate()
        v_rhs = (- u + alpha + D*dxdxv).evaluate()   
        v_rhs.change_scales(1)
        v_rhs_list.append(v_rhs['g'].copy())

        while solver.proceed and solver.sim_time < t_stop - 0.5*timestep:
            solver.step(timestep)

            if solver.iteration % log_every == 0:
                logger.info(
                    "Iteration=%i, Time=%e, dt=%e, norm(u)=%e"
                    % (solver.iteration, solver.sim_time, timestep, np.linalg.norm(u['g']))
                )

            if solver.iteration % save_every == 0:
                u.change_scales(1)
                v.change_scales(1)

                u_list.append(u['g'].copy())
                v_list.append(v['g'].copy())
                t_list.append(solver.sim_time)
                
                dxu = d3.Differentiate(u, xcoord).evaluate()
                dxdxu = d3.Differentiate(dxu, xcoord).evaluate()
                u_rhs = ((v - u*u - u*u*u)/eps + D*dxdxu).evaluate()

                dxv = d3.Differentiate(v, xcoord).evaluate()
                dxdxv = d3.Differentiate(dxv, xcoord).evaluate()
                v_rhs = (- u + alpha + D*dxdxv).evaluate()
                
                u_rhs.change_scales(1)
                v_rhs.change_scales(1)

                u_rhs_list.append(u_rhs['g'].copy())
                v_rhs_list.append(v_rhs['g'].copy())

        t_arr = np.array(t_list)     
        u_snap = np.array(u_list)    
        v_snap = np.array(v_list)    
        u_rhs_snap = np.array(u_rhs_list)   
        v_rhs_snap = np.array(v_rhs_list)   

        np.savez(out_npz, 
            x=x.copy(), 
            t=t_arr, 
            u=u_snap, 
            v=v_snap, 
            u_rhs=u_rhs_snap, 
            v_rhs=v_rhs_snap, 
            dt=float(timestep * save_every)
        )
        logger.info(f"Saved full trajectory to {out_npz}")

    if POD:
        # ---- load data ----
        data = np.load(out_npz)

        x = data["x"].copy()
        U = data["u"].T          # (Nx, N)
        V = data["v"].T          # (Nx, N)
        U_rhs = data["u_rhs"].T  # (Nx, N)
        V_rhs = data["v_rhs"].T  # (Nx, N)
        t = data["t"]
        dt = data["dt"]

        Nx, N = U.shape

        print(f"N_train: {N_train}, N_test: {N_test}")
        train_idx = np.arange(N_train)
        test_idx = np.arange(N_train, N_train + N_test)

        # ---- split train / test ----
        U_train = U[:, train_idx]
        V_train = V[:, train_idx]
        U_test  = U[:, test_idx]
        V_test  = V[:, test_idx]

        U_rhs_train = U_rhs[:, train_idx]
        V_rhs_train = V_rhs[:, train_idx]
        U_rhs_test  = U_rhs[:, test_idx]
        V_rhs_test  = V_rhs[:, test_idx]

        t_train = t[train_idx]
        t_test  = t[test_idx]

        # ---- quadrature weights ----
        W = compute_quadrature_weights(dist, xbasis, Nx)   # (Nx,)
        # print("sum(W) =", np.sum(W))   # should be about 2 on [-1,1]

        # JOINT POD on stacked state Y = [U; V] in R^{2Nx x N}
        Y_train = np.vstack([U_train, V_train])            # (2Nx, N_train)
        Y_test  = np.vstack([U_test,  V_test])             # (2Nx, N_test)

        # stacked RHS snapshots = time derivative of stacked state
        Y_rhs_train = np.vstack([U_rhs_train, V_rhs_train])   # (2Nx, N_train)
        Y_rhs_test  = np.vstack([U_rhs_test,  V_rhs_test])    # (2Nx, N_test)

        # quadrature weights for joint state: repeat W for U and V components
        W_joint = np.concatenate([W, W])                  # (2Nx,)

        for n_pod_mode in n_pod_modes:
            # ---- weighted POD of stacked snapshots on TRAIN only ----
            Phi, Y_mean, xi_train, sing_vals = weighted_pod(
                Y_train, W_joint, r=n_pod_mode, subtract_mean=True
            )

            if plot_singular_values:
                plt.figure()
                plt.semilogy(np.arange(n_pod_mode), sing_vals[:n_pod_mode], 'o-')
                plt.xlabel("POD modes")
                plt.ylabel("singular value")
                plt.title("Singular values from weighted POD")
                plt.tight_layout()
                # plt.show()
                plt.savefig(f"{plot_dir}/pod_singular_values_Nx={Nx}_d={n_pod_mode}_N={N_total}.png", dpi=300)
                plt.close()
            
            # ensure mean has shape (2Nx, 1)
            Y_mean = Y_mean.reshape(2 * Nx, 1)

            # ---- reduced coefficients on TEST ----
            xi_test = Phi.T @ (W_joint.reshape(2 * Nx, 1) * (Y_test - Y_mean))

            # ============================================================
            # TRUE reduced dynamics: xi_dot = Phi^T W_joint Y_rhs
            # ============================================================
            xi_dot_test  = Phi.T @ (W_joint.reshape(2 * Nx, 1) * Y_rhs_test)

            # ============================================================
            # Reconstruction of test data from true POD coefficients
            # (useful baseline before comparing with learned xi_hat)
            # ============================================================
            plot_reconstruction = False
            if plot_reconstruction:
                Y_rec_test  = Y_mean + Phi @ xi_test             # (2Nx, N_test)
                u_hat_test = Y_rec_test[:Nx, :]
                v_hat_test = Y_rec_test[Nx:, :]

                # plot reconstruction error for train and test 
                plt.figure(figsize=(10, 4))
                plt.subplot(1, 2, 1)
                plt.pcolor(t_test, x, np.abs(U_test - u_hat_test), label="test U")
                plt.colorbar()
                plt.xlabel("t")
                plt.ylabel("x")
                plt.title("abs. error on test set for u")

                plt.subplot(1, 2, 2)
                plt.pcolor(t_test, x, np.abs(V_test - v_hat_test), label="test V")
                plt.colorbar()
                plt.xlabel("t")
                plt.ylabel("x")
                plt.title("abs. error on test set for v")
                plt.tight_layout()
                plt.savefig(f"{plot_dir}/reconstruction_error_Nx={Nx}_d={n_pod_mode}_N={N_total}.png", dpi=300)
                plt.show()



            # # Save 
            out_h5_train = f"{data_dir}/PODcoeffs_d={n_pod_mode}_N={N_total}_train.hdf5"
            out_h5_test = f"{data_dir}/PODcoeffs_d={n_pod_mode}_N={N_total}_test.hdf5"


            
            with h5py.File(out_h5_train, "w") as f:
                f.create_dataset("x", data=xi_train.T)   # reduced state
                f.create_dataset("dt", data=float(dt))

            # Save test file
            with h5py.File(out_h5_test, "w") as f:
                f.create_dataset("x", data=xi_test.T)         # reduced state
                f.create_dataset("f", data=xi_dot_test.T)     # reduced dynamics

                # optional extra fields
                f.create_dataset("x_grid", data=x.copy())     # renamed spatial variable
                f.create_dataset("t", data=t.copy())
                f.create_dataset("Phi", data=Phi.copy())
                f.create_dataset("Y_mean", data=Y_mean.copy())
                f.create_dataset("u_test", data=U_test.copy())
                f.create_dataset("v_test", data=V_test.copy())

            print(f"Saved POD coefficients to {out_h5_train} and {out_h5_test}")
