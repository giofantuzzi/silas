import numpy as np
import dedalus.public as d3

# -----------------------
# Chebyshev T_n fields and coefficient computation
# -----------------------
def make_Tn_fields(dist, xbasis, x, orders):
    xi = x  
    Tn_fields = []
    for n in orders:
        Tn = dist.Field(name=f"T{n}", bases=xbasis)
        Tn['g'] = np.cos(n * np.arccos(np.clip(xi, -1.0, 1.0)))  # T_n(xi)
        Tn_fields.append(Tn)
    return Tn_fields


def cheb_coeffs(u, Tn_fields):
    # a_n(t) = \int T_n(x) u(x,t) dx
    a = []
    for Tn in Tn_fields:
        val = d3.Integrate(Tn*u, 'x').evaluate()['g']
        a.append(np.asarray(val).item())  # explicit scalar extraction
    return np.array(a)


# -----------------------
# Running trapezoid means
# -----------------------
class RunningMeans:
    def __init__(self, m):
        self.m = m
        self.t_prev = None
        self.a_prev = None
        self.I = np.zeros(m)  # integral of a_n
        self.t = []
        self.a = []
        self.mean = []

    def update(self, t, a_vec):
        if self.t_prev is None:
            self.t_prev = t
            self.a_prev = a_vec.copy()
            self.t.append(t)
            self.a.append(a_vec.copy())
            self.mean.append(np.zeros_like(a_vec))
            return

        dtloc = t - self.t_prev
        self.I += 0.5 * (a_vec + self.a_prev) * dtloc

        self.t_prev = t
        self.a_prev = a_vec.copy()

        self.t.append(t)
        self.a.append(a_vec.copy())
        self.mean.append(self.I / t)

    def arrays(self):
        return np.array(self.t), np.array(self.a), np.array(self.mean)



def compute_quadrature_weights(dist, xbasis, Nx):
    tmp = dist.Field(bases=xbasis)
    tmp.change_scales(1)
    x = dist.local_grid(xbasis, scale=1)

    basis_eval = np.zeros((Nx, Nx))
    modal_int = np.zeros(Nx)

    for n in range(Nx):
        tmp['c'] = 0
        tmp['c'][n] = 1
        tmp.change_scales(1)
        basis_eval[:, n] = tmp['g'].copy()
        modal_int[n] = d3.Integrate(tmp, 'x').evaluate()['g'].item()

    w = np.linalg.solve(basis_eval.T, modal_int)
    return w



def weighted_pod(U: np.ndarray, W: np.ndarray, r: int, subtract_mean: bool = True):
    """
    U: (Nx, K) snapshot matrix for u
    W: (Nx,) quadrature weights
    r: number of POD modes

    Returns:
        Phi: (Nx, r) L2-orthonormal POD modes, Phi^T W Phi = I
        u_mean: (Nx, 1) mean used (zeros if subtract_mean=False)
        coeffs: (r, K) POD coefficients: Phi^T W (U - mean)
        S: singular values from SVD of W^{1/2} (U - mean)
    """
    Nx, K = U.shape

    W = W.reshape(-1)
    assert W.shape == (Nx,)

    if subtract_mean:
        u_mean = np.mean(U, axis=1, keepdims=True)
        Uc = U - u_mean
    else:
        u_mean = np.zeros((Nx, 1), dtype=U.dtype)
        Uc = U

    wsqrt = np.sqrt(W).reshape(Nx, 1)
    winvsqrt = (1.0 / np.sqrt(W)).reshape(Nx, 1)

    # U_tilde = W^{1/2} Uc
    U_tilde = wsqrt * Uc

    # SVD 
    Uu, S, Vh = np.linalg.svd(U_tilde, full_matrices=False)

    # Phi = W^{-1/2} U_tilde
    Phi = winvsqrt * Uu[:, :r]

    # Coeffs: \Xi = \Phi^T W Uc
    coeffs = Phi.T @ (W.reshape(Nx, 1) * Uc)

    return Phi, u_mean, coeffs, S