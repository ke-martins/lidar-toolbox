# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars 🌊

This is part of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidars. This sub-directory **lidBathy** contains a complete workflow example for the implementation of the non-linear depth inversion algorithm of Martins et al. (2023) on field lidar data. It relies on most libraries within **lidar-toolbox** to grid the raw multibeam lidar data (**lidar-gridding**), and perform spectral, cross-spectral and bispectral analyses of gappy free surface elevation timeseries (**gappy-series-preprocessing** and **spectral-analysis**). This example accompanies the manuscript describing this lidar implementation of *lidBathy* and submitted to *CENG* (Martins et al., submitted). A list of other relevant publications that used this toolbox is given at the bottom of this page.  

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
spectral-analysis v1: first release of the library.

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

# References

 - Martins, K., Bonneton, P., & Michallet, H. (2021). Dispersive characteristics of non-linear waves propagating and breaking over a mildly sloping laboratory beach. *Coastal Engineering* <strong>167</strong>, 103917. https://doi.org/10.1016/j.coastaleng.2021.103917
 
 - Martins, K., Bonneton, P., Lannes, D., & Michallet, H. (2021). Relation between orbital velocities, pressure, and surface elevation in nonlinear nearshore water waves. *Journal of Physical Oceanography* <strong>51</strong>(11), 3539-3556. https://doi.org/10.1175/JPO-D-21-0061.1
 
 - Martins, K., Bonneton, P., de Viron, O., Turner, I. L., Harley, M. D., & Splinter, K. (2023). New Perspectives for Nonlinear Depth‐Inversion of the Nearshore Using Boussinesq Theory. *Geophysical Research Letters* <strong>50</strong>(2), e2022GL100498. https://doi.org/10.1029/2022GL100498
 
 - Martins, K., Brodie, K. L., Fiedler, J. W., O'Dea, A. M., Spore, N. J., Grenzeback, R. L., Dickhudt, P. J., Bak, S. de Viron, O., Bonneton, P. Seamless nearshore topo-bathymetry reconstruction from lidar scanners: a Proof-of-Concept based on a dedicated field experiment at Duck, NC. *submitted to Coastal Engineering*.
 
 - Sous, D., Martins, K., Tissier, M., Bouchette, F., & Meulé, S. (2023). Spectral wave dissipation over a roughness‐varying barrier reef. *Geophysical Research Letters* <strong>50</strong>(5), e2022GL102104. https://doi.org/10.1029/2022GL102104
 
 
 
 
 
 
