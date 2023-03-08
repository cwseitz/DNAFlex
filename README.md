# DNA single molecule localization microscopy (DNA-SMLM)

Abbreviated as just SMLM, is a library built on Python and C for SMLM experiment design, optimization, and analysis

## Current modules

A toolkit for simulations and analysis related to DNA SMLM in Python. The package will include: 

1. **localize** Basic particle detection and tracking functions for static and dynamic single molecule data
2. **psf** Gaussian point spread function models in 3D, sCMOS noise models, and information theoretic localization techniques
3. **switch** Functions for performing Monte Carlo simulations of photoswitching dynamics of fluorophores
4. **plot** A variety of functions for generating common SMLM plots and making animations

## Modules in progress

5. **md** Coarse-grained molecular dynamics simulations for nucleosome ball and stick type models
6. **stats** Statistical software e.g., Markov Chain Monte Carlo samplers for performing Bayesian inference on SMLM-related models, maximum likelihood estimators
7. **torch** A place for any PyTorch deep models related to SMLM

Certain modules with a subfolder **_MODULE** contain backend C code for optimization.
