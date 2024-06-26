# extrapolated-regularization
Code for https://arxiv.org/abs/2309.14169

For questions/comments contact Svetlana Tlupova (svetlana.tlupova@farmingdale.edu) 

Reference:
J. T. Beale and S. Tlupova. Extrapolated regularization of nearly singular integrals on surfaces. arXiv; Cornell University Library, 2024. https://arxiv.org/abs/2309.14169

This code was developed as part of work supported by the National Science Foundation grant DMS-2012371.

## Testing instructions

### Test: Harmonic potentials on a molecular surface (Fig. 9 in referenced paper, Table 10 in original version on arxiv):

1.	Code subdirectory: `extrap_harmonic_molecule`
2.	Run Matlab code

### Test: Stokeslet and stresslet integrals on a molecular surface (Fig. 13 in referenced paper, Table 16 in original version on arxiv):

1.	Code subdirectory: `extrap_stokes_molecule`
2.	In utilities.h, set the following parameters:
    * grid size (h=1/N_gridlines): N_gridlines = 32, 64, 128, 256
    * type of regularization in constant “del_h”:
        * del_h = 1:  delta=c*h
        * del_h = 2 (anything not equal 1):  delta=c*h^(⅘)
3.  Compile using `make`
4.  Run using `./main.out`
