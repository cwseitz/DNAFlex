# DNAFlex

Abbreviated as just DNAFlex, this is a library built on Python and C for single molecule localization microscopy (SMLM) experiment design, optimization, and analysis associated with DNA structure at the scale of rigid nucleosomes

## Current modules

A toolkit for simulations and analysis related to DNA SMLM in Python. The package will include: 

1. **localize** Basic particle detection and tracking functions for static and dynamic single molecule data
2. **psf** Gaussian point spread function models in 3D, sCMOS noise models, and information theoretic localization techniques
3. **ssa** Functions for performing Monte Carlo simulations of photoswitching dynamics of fluorophores
4. **plot** A variety of functions for generating common SMLM plots and making animations

## Modules in progress

5. **flex** DNA Flex: Monte Carlo simulations of flexible linker DNA and rigid nucleosomes
6. **stats** Statistical software e.g., Markov Chain Monte Carlo samplers for performing Bayesian inference on SMLM-related models, maximum likelihood estimators
7. **torch** A place for any PyTorch deep models related to SMLM

Certain modules with a subfolder **_MODULE** contain backend C code for optimization.

## How to run DNAFlex C code

```
sh run.sh nucl_pos.dat link_seq.dat 10 helpars.dat /home/cwseitz/Desktop/test
```
