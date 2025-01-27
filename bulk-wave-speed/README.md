# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars

This is part of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidars. This sub-directory **bulk-wave-speed** contains only one function that computes the mean wave propagation speed from two free surface elevation signals.

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
bulk-wave-speed v1: first release of the library, no example for now.

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

## List of functions

<details>
  <summary>📄 fun_compute_c_from_xcorr.m</summary>  
  <br>  

  **Description**:  
  Function computing 'bulk' wave propagation speed from cross-correlation of adjacent surface elevation timeseries. We use spline interpolation to reach sub-timestep resolution. For now, this assumes no NaNs in the provided timeseries so these have to be dealt with elsewhere.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `zeta_2`  | double | surface elevation timeseries at position 2 [m]                  |
  | `zeta_1`  | double | surface elevation timeseries at position 1 [m]                  |
  | `sf`  | int | sampling frequency [Hz]                 |
  | `dx`  | double | distance between two locations where zeta is provided [m]               |
  | `maxlag`  | double | max lag authorised for cross-correlation analysis [s]                |

  **Outputs**:  
  &nbsp;&nbsp;`c`, mean wave speed [m/s]  
  &nbsp;&nbsp;`lag`, corresponding time lag [s]  

</details>


 
 
