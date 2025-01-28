# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars 🌊

Welcome to this short description of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidar scanners. The two main features are the libraries for gridding multibeam lidar (**lidar-gridding**) and performing spectral, cross-spectral and bispectral analysis (**spectral-analysis**). In essence, the **spectral-analysis** library is an adaptation of my **bispectral-analysis** library to gappy data collected by multibeam lidar scanners (https://github.com/ke-martins/bispectral-analysis). This adaptation was motivated while developing the nonlinear, lidar-based, nearshore depth inversion algorithm and working on the manuscript submitted to *CENG* (Martins et al., submitted).  

At the moment, there are 4 inter-dependent sub-libraries, each dedicated to specific functions described below, but this is expected to grow in the future, and include libraries for data acquisition for instance. In the examples folder, there several example scripts to grid multibeam lidar data and perform spectral analysis on such gridded data. Additionally, a complete workflow example of the implementation of the Boussinesq-based neashore depth inversion method of Martins et al. (2023) to field lidar datasets is presented. This example is released as a accompanying code for the manuscript submitted to *CENG* (Martins et al., see references list).  

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
lidar-toolbox v1: first release of the library. Fully tested and examples provided.

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

## List of directories/sub libraries

<details>
  <summary>📂 <strong>lidar-gridding</strong> </summary>  
  <br>  

  **Description**:  
  Library containing functions that grid (spatial and temporal interpolations) data collected by single- and multibeam lidar scanners. It is based on the built-in function `scatteredInterpolant` (Delaunay triangulation) and makes the most of 4D (t,x,y,z) point clouds to interpolate raw data on grids and respect some quality criteria to prevent gaps over-filling.

  **List of functions**:  

  fun_gridded_lidar_diagnostics.m  
  fun_singlebeam_lidar_gridding.m  
  fun_multibeam_lidar_gridding.m

</details>

<details>
  <summary>📂 <strong>gappy-series-preprocessing</strong> </summary>  
  <br>  

  **Description**:  
  Functions for pre-processing gappy data series so that spectral, cross-spectral and bispectral analyses can be applied to them. Series are essentially reorganised by blocks, and quality is controlled through the number of NaNs allowed per block.

  **List of functions**:  

  fun_count_pNaNs.m  
  fun_interp_series.m  
  fun_prep_gappy_series_by_block.m
  fun_prep_gappy_series_by_block_xs.m  

</details>

<details>
  <summary>📂 <strong>spectral-analysis</strong> </summary>  
  <br>  

  **Description**:  
  Functions needed to perform spectral, cross-spectral and bispectral analyses on gappy free surface elevation timeseries of ocean waves measured with lidars. Timeseries should be pre-organised in matrices with the library **gappy-series-preprocessing**. For users interested in bispectral products, it also directly contains relevant functions for a range of nearshore applications (wave dispersive properties, non-linear energy transfers between triads etc).

  **List of functions**:  

  fun_compute_spectrum_mat.m  
  fun_compute_cross_spectrum_mat.m  
  fun_compute_bispectrum_mat.m
  fun_compute_krms.m  
  fun_compute_krms_terms.m  
  fun_compute_Snl.m  
  fun_compute_edof.m  

</details>

<details>
  <summary>📂 <strong>bulk-wave-speed</strong> </summary>  
  <br>  

  **Description**:  
  Function to compute sub-resolution lag between two timeseries and corresponding bulk celerity.

  **List of functions**:  

  fun_compute_c_from_xcorr.m  

</details>

# References
 
 - Martins, K., Bonneton, P., de Viron, O., Turner, I. L., Harley, M. D., & Splinter, K. (2023). New Perspectives for Nonlinear Depth‐Inversion of the Nearshore Using Boussinesq Theory. *Geophysical Research Letters* <strong>50</strong>(2), e2022GL100498. https://doi.org/10.1029/2022GL100498
 
 - Martins, K., Brodie, K. L., Fiedler, J. W., O'Dea, A. M., Spore, N. J., Grenzeback, R. L., Dickhudt, P. J., Bak, S. de Viron, O., Bonneton, P. Seamless nearshore topo-bathymetry reconstruction from lidar scanners: a Proof-of-Concept based on a dedicated field experiment at Duck, NC. *submitted to Coastal Engineering*.  
 
 
 
 
 
 
