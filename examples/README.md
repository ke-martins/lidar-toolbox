# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars 🌊

This is part of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidars. This sub-directory contains several examples to handle lidar data, including a complete workflow example for the implementation of the non-linear depth inversion algorithm of Martins et al. (2023) on field lidar data. It relies on most libraries within **lidar-toolbox** to grid the raw multibeam lidar data (**lidar-gridding**), and perform spectral, cross-spectral and bispectral analyses of gappy free surface elevation timeseries (**gappy-series-preprocessing** and **spectral-analysis**). This example accompanies the manuscript describing this lidar implementation of *lidBathy* and submitted to *CENG* (Martins et al., submitted). A list of other relevant publications that used this toolbox is given at the bottom of this page.  

## 📥 Please download the pre-processed file before running example

### **<sub><sup>1️⃣</sup></sub> Manual Download**
Follow the link below, download the required file and place it into data/ directory:  
- [Download 20220912_191123_flight_2_x=225.mat](https://drive.google.com/uc?id=13Qk01nyErT1LT-RtmG7dNQtoTl9gqpHI) 📄

---

### **<sub><sup>2️⃣</sup></sub> Directly through the MATLAB code provided and the `gdown`**
The MATLAB code will install the gdown command if needed, and download directly the file (Unix systems only). Please use option 1 is you do not wish to install the gdown command.  
#### **Using `curl`**
```bash
pip install gdown
gdown https://drive.google.com/uc?id=13Qk01nyErT1LT-RtmG7dNQtoTl9gqpHI

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
