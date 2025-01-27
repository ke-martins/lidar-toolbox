# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars 🌊

Welcome to this short description of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidar scanners. The two main features are the libraries for gridding multibeam lidar (**lidar-gridding**) and performing spectral, cross-spectral and bispectral analysis (**spectral-analysis**). In essence, the **spectral-analysis** library is an adaptation of my **bispectral-analysis** library to lidar data (https://github.com/ke-martins/bispectral-analysis). This adaptation was motivated while developing the nonlinear, lidar-based, nearshore depth inversion algorithm and working on the manuscript submitted to *CENG* (Martins et al., submitted). A list of other relevant publications that used this toolbox is given at the bottom of this page. Although this toolbox was intended for surface elevation datasets, it can be used for other signals bearing in mind the provided units would be wrong. At the moment, there are 4 inter-dependent sub-libraries, each dedicated to specific functions described below, but this is expected to grow in the future, and include libraries for data acquisition for instance.  

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
lidar-toolbox v1: first release of the library.

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

## List of directories/sub libraries

<details>
  <summary>📂 <strong>lidar-gridding<strong> </summary>  
  <br>  

  **Description**:  
  Library containing functions that grid (spatial and temporal interpolations) data collected by single- and multibeam lidar scanners. It is based on the function built function `scatteredInterpolant` (Delaunay triangulation) and makes the most of 4D (t,x,y,z) point clouds to interpolate raw data on grids and respect some quality criteria.

  **List of functions**:  

  fun_gridded_lidar_diagnostics.m  
  fun_singlebeam_lidar_gridding.m  
  fun_multibeam_lidar_gridding.m

</details>

<details>
  <summary>📂 <strong>gappy-series-preprocessing<strong> </summary>  
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
  <summary>📂 <strong>spectral-analysis<strong> </summary>  
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
  <summary>📂 <strong>bulk-wave-speed<strong> </summary>  
  <br>  

  **Description**:  
  Function to compute sub-resolution lag between two timeseries and corresponding bulk celerity.

  **List of functions**:  

  fun_compute_c_from_xcorr.m  

</details>

# References

 - Martins, K., Bonneton, P., & Michallet, H. (2021). Dispersive characteristics of non-linear waves propagating and breaking over a mildly sloping laboratory beach. *Coastal Engineering* <strong>167</strong>, 103917. https://doi.org/10.1016/j.coastaleng.2021.103917
 
 - Martins, K., Bonneton, P., Lannes, D., & Michallet, H. (2021). Relation between orbital velocities, pressure, and surface elevation in nonlinear nearshore water waves.* Journal of Physical Oceanography* <strong>51</strong>(11), 3539-3556. https://doi.org/10.1175/JPO-D-21-0061.1
 
 - Martins, K., Bonneton, P., de Viron, O., Turner, I. L., Harley, M. D., & Splinter, K. (2023). New Perspectives for Nonlinear Depth‐Inversion of the Nearshore Using Boussinesq Theory. *Geophysical Research Letters* <strong>50</strong>(2), e2022GL100498. https://doi.org/10.1029/2022GL100498
 
 - Martins, K., Brodie, K. L., Fiedler, J. W., O'Dea, A. M., Spore, N. J., Grenzeback, R. L., Dickhudt, P. J., Bak, S. de Viron, O., Bonneton, P. Seamless nearshore topo-bathymetry reconstruction from lidar scanners: a Proof-of-Concept based on a dedicated field experiment at Duck, NC. *submitted to Coastal Engineering*.
 
 - Sous, D., Martins, K., Tissier, M., Bouchette, F., & Meulé, S. (2023). Spectral wave dissipation over a roughness‐varying barrier reef. *Geophysical Research Letters* <strong>50</strong>(5), e2022GL102104. https://doi.org/10.1029/2022GL102104
 
 
 
 
 
 
