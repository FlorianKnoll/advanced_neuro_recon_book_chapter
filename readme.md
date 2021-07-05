-----------------------------------------------
Matlab code for the didactic experiments for the book chapter "Advanced Reconstruction Methods for Fast MRI" for Advanced Neuro MR Techniques and Applications
June 2020
Florian Knoll (florian.knoll@nyumc.org)
-----------------------------------------------

This package contains source code for most of the didactic examples in the image reconstruction chapter of the textbook Advanced Neuro MR Techniques and Applications. The following scripts recreate the reconstruction examples for the chapters partial Fourier imaging, parallel imaging with SENSE and GRAPPA, and Compressed Sensing. The training code for the machine learning chapter is available at: https://github.com/VLOGroup/mri-variationalnetwork, and the test set reconstruction results are provided for reference.

The data and pre-estimated coil sensitivity maps can be obtained via globus at the following location:
Link: https://app.globus.org/file-manager?origin_id=15c7de28-a76b-11e9-821c-02b7a92d8e58&origin_path=%2F
Subfolder: advanced_neuro_mr_techniques_applications/

-----------------------------------------------
pf_pocs_reconstruction_example.m
-----------------------------------------------
The script simulates a 5/8 half Fourier acquisition and shows reconstructions from simple zero filling, brute force conjugate symmetry reconstruction as well as an iterative phase constrained reconstruction using the projection on convex sets (POCS) algorithm [1]. Please note that the original dataset is from a TSE acquisition and does not show severe residual phase. Therefore, I modulated it with an additionally synthetic sinusoidal phase to enhance the effects of phase errors, and added Gaussian noise to the k-space for didactic reasons. This is part of the script, so if you want to see the partial Fourier reconstruction without this additional corruption, you can just comment this section out. 

----------------------------------------------
SENSE_reconstruction_example.m
-----------------------------------------------
Iterative (CG) SENSE reconstruction based on [2]. Regular equidistant sampling masks are provided for R=4 and 26 reference lines at the center of k-space and R=6 and 16 lines. Coil sensitivities were pre-estimated using ESPIRIT [3]. The core CG-SENSE algorithm is implemented in the function pmri_cgsense.m in the subfolder sense.

-------------------------------------------------
GRAPPA_reconstruction_example.m
-------------------------------------------------
GRAPPA reconstruction based on [4]. Regular equidistant sampling masks are provided for R=4 and 26 reference lines at the center of k-space and R=6 and 16 lines. The  implementation for this didactic example is based on the GRAPPA implementation that is provided as part of Miki Lustig's code package for SPIRiT [5], which is available online from https://people.eecs.berkeley.edu/~mlustig/softwaregrappa/SPIRiT_v0.1.tar.gz
The following functions from this package need to be added to this folder: calibrate.m, corrMatrix.m, getCalibSize.m, GRAPPA.m

-------------------------------------------
TGV_reconstruction_example.m
-------------------------------------------
Combined parallel imaging and compressed sensing reconstruction with a Total Generalized Variation (TGV) constraint [6]. Pseudo-random sampling masks are provided for R=4 and 26 reference lines at the center of k-space and R=6 and 16 lines.The implementation follows the code that was provided for that paper, it was modified for Cartesian acquisitions as opposed to radial acquisitions in the original code of the MRM paper. The core optimization is a first order primal dual algorithm [7], implemented in the function tgv2_l2_2D_pd.m.

----------------
References
----------------
[1] Cuppen J, van Est A. Reducing MR imaging time by one-sided reconstruction. Magnetic Resonance Imaging. 5:526-527 (1987).
[2] Pruessmann KP, Weiger M, Boernert P, Boesiger P. Advances in sensitivity encoding with arbitrary k-space trajectories. Magnetic Resonance in Medicine 46: 638-651 (2001).
[3] Uecker M, Lai P, Murphy MJ, Virtue P, Elad M, Pauly JM, Vasanawala SS, Lustig M. ESPIRiT - An eigenvalue approach to autocalibrating parallel MRI: Where SENSE meetsGRAPPA. Magnetic Resonance in Medicine 71:, 990–1001 (2014).
[4] Griswold MA, Jakob PM, Heidemann RM, Nittka M, Jellus V, Wang J, Kiefer B, Haase A. Generalized Autocalibrating Partially Parallel Acquisitions (GRAPPA). Magnetic Resonancein Medicine 47: 1202-1210 (2002).
[5] Lustig M, Pauly J. SPIRiT: Iterative Self-Consistent Parallel Imaging Reconstruction from Arbitrary k-Space Sampling. Magnetic Resonance in Medicine. 64:457-471 (2010). 
[6] Knoll F, Bredies K, Pock T, Stollberger, R. Second Order Total Generalized Variation (TGV) for MRI. Magnetic Resonance in Medicine. 65: 480-491 (2011).
[7] Chambolle A, Pock T. A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging. J Math Imaging Vis 40: 120–145 (2011).
