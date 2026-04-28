# Data-driven discovery of bounded polynomial ODEs

This repository contains code and data used to produce results in the following paper:

[paper-here]

Most of the code is in MATLAB, with some Python scripts to generate data and post-process results.

The code is provided "as-is" and without any guarantees. See the license file for more information.

**Authors:** Albert Alcalde, Giovanni Fantuzzi (FAU Erlangen-Nuernberg)

## Dependencies
To run the code, you must first install the following MATLAB packages:
* YALMIP
* MOSEK (or another SDP solver compatible with YALMIP)
* ChebFun

## Cite us
If you find this code useful in your work and publish your results, please cite the arXiv preprint above. An sample .bib code follows.

```
    @article{FantuzziAlcalde2026,
        AUTHOR = {Alcalde, Albert and Fantuzzi, Giovanni},
        TITLE = {Data-driven discovery of polynomial ODEs with provably bounded solutions},
        YEAR = {2026},
        HOWPUBLISHED = {arXiv:????.????},
        URL = {https://arxiv.org/abs/????.????}
    }
```
