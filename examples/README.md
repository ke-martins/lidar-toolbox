# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars 🌊

This is part of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidars. This sub-directory contains several examples to handle lidar data, including a complete workflow example for the implementation of the non-linear depth inversion algorithm of Martins et al. (2023) on field lidar data. It relies on most libraries within **lidar-toolbox** to grid the raw multibeam lidar data (**lidar-gridding**), and perform spectral, cross-spectral and bispectral analyses of gappy free surface elevation timeseries (**gappy-series-preprocessing** and **spectral-analysis**). This example accompanies the manuscript describing this lidar implementation of *lidBathy* and submitted to *CENG* (Martins et al., submitted). A list of other relevant publications that used this toolbox is given at the bottom of this page.  

### Prerequisite: pre-processed lidar data file for SIO flight #2 of BELS2022 experiments  
Please make sure you download the pre-processed lidar data file before running example. It is pretty heavy, and hosted elsewhere as github do not easily accommodate for files this big. Two options are proposed below:  
&nbsp;&nbsp;<sub><sup>1️⃣</sup></sub> Manually download the file at the following link; and place it into `data/` directory: [20220912_191123_flight_2_x=225.mat](https://drive.google.com/uc?id=13Qk01nyErT1LT-RtmG7dNQtoTl9gqpHI) 📄  
&nbsp;&nbsp;<sub><sup>2️⃣</sup></sub> Directly through the MATLAB code provided and the `gdown` command  
Please use option 1 is you do not wish to install the `gdown` command, though it is a super light utility. The MATLAB code (multibeam_lidar_gridding.m or lidBathy_workflow_example.m) will install the `gdown` utility if needed, and download the file through the following commands (Unix systems only for now):    
```bash
pip install gdown
gdown https://drive.google.com/uc?id=13Qk01nyErT1LT-RtmG7dNQtoTl9gqpHI
mv 20220912_191123_flight_2_x=225.mat data/
```

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
spectral-analysis v1: first release of the library.

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

# References
 
 - Martins, K., Bonneton, P., de Viron, O., Turner, I. L., Harley, M. D., & Splinter, K. (2023). New Perspectives for Nonlinear Depth‐Inversion of the Nearshore Using Boussinesq Theory. *Geophysical Research Letters* <strong>50</strong>(2), e2022GL100498. https://doi.org/10.1029/2022GL100498
 
 - Martins, K., Brodie, K. L., Fiedler, J. W., O'Dea, A. M., Spore, N. J., Grenzeback, R. L., Dickhudt, P. J., Bak, S. de Viron, O., Bonneton, P. Seamless nearshore topo-bathymetry reconstruction from lidar scanners: a Proof-of-Concept based on a dedicated field experiment at Duck, NC. *submitted to Coastal Engineering*.  
