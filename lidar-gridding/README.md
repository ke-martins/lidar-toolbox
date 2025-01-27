# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars

This is part of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidars. This sub-directory **lidar-gridding** contains functions that grid (spatial and temporal interpolations) data collected by single- or multibeam lidar scanners. It is based on the function built function `scatteredInterpolant` (Delaunay triangulation) and makes the most of 4D (t,x,y,z) point clouds to interpolate raw data on grids and respect some quality criteria.

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
lidar-gridding v1: first release of the library; example to be added separately, while reseasing entire lidBathy workflow.

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

## List of functions

<details>
  <summary>📄 fun_gridded_lidar_diagnostics.m</summary>  
  <br>  

  **Description**:  
  Computing basic statistics of gridded lidar data.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `sf`    | double | sampling frequency [Hz]                    |
  | `x`    | double | cross-shore grid [m]                |
  | `data`    | double | gridded data (dimensions: (t,x)), typically the field 'z' output by fun_multibeam_lidar_gridding                   |
  | `bplot`    | int | optional, input for plotting: 1 for yes, 0 for no (default)              |

  **Outputs**:  
  &nbsp;&nbsp;Returns `s`, a self-explanatory data structure containing some basic info on gridded lidar `data`.

</details>

<details>
  <summary>📄 fun_singlebeam_lidar_gridding.m</summary>  
  <br>  

  **Description**:  
  Function gridding (both time and space) the raw data from singlebeam lidar systems. This was initially written for processing the data collected at Duck with UAV-mounted lidar systems. It assumes minimal pre-processing (and filtering) for the input data, with the idea to make this function as generic as possible. At minima, it needs a time (array of size Np x 1) and xyz (array of size Np x 3) to work properly.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `time`    | double | time interpolation grid |  
  | `x_grid`  | double | x interpolation grid [m]   |  
  | `raw_data` | struct | data structure containing de-noised raw point cloud; minimum data fields: time (Np x 1); xyz (Np x 3) |  
  | `t_win`  | double | time window [s] within which data is used for the time interpolation; e.g., if t_win = 0.15 s, then for each time ti of the grid, the algorithm will use data within 0.15 s both in past (ti-0.15) and future (ti+0.15). |  
  | `x_win`  | double | cross-shore space window [m] within which data needs to be present for keeping interpolated value (NaN is used otherwise) |  

  **Outputs**:  
  &nbsp;&nbsp;`grid_data`, a self-explanatory data structure containing the gridded lidar data.  
  
  **Comments on parameters `t_win` and `x_win`**:  
  In some way, `t_win` and `x_win` are used to select the interpolation data points and indirectly define regions outside of which, data is not interpolated. So increasing those two parameters will fill in more gaps, but potentially increase the number of irrealistic points. This is slightly less true for `t_win`, but keep in mind that increasing `t_win` also slows down the interpolation process since more data is being used to create the interpolant.  

</details>

<details>
  <summary>📄 fun_multibeam_lidar_gridding.m</summary>  
  <br>  

  **Description**:  
  Function gridding (both time and space) the raw data from singlebeam lidar systems. This was initially written for processing the data collected at Duck with UAV-mounted lidar systems. It assumes minimal pre-processing (and filtering) for the input data, with the idea to make this function as generic as possible. At minima, it needs a time (array of size Np x 1) and xyz (array of size Np x 3) to work properly.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `time`    | double | time interpolation grid |  
  | `x_grid`  | double | x interpolation grid [m]   |  
  | `y_grid`  | double | y interpolation grid [m]   |  
  | `raw_data` | struct | data structure containing de-noised raw point cloud; minimum data fields: time (Np x 1); xyz (Np x 3) |  
  | `t_win`  | double | time window [s] within which data is used for the time interpolation; e.g., if t_win = 0.15 s, then for each time ti of the grid, the algorithm will use data within 0.15 s both in past (ti-0.15) and future (ti+0.15). |  
  | `x_win`  | double | cross-shore space window [m] within which data needs to be present for keeping interpolated value (NaN is used otherwise) |  

  **Outputs**:  
  &nbsp;&nbsp;`grid_data`, a self-explanatory data structure containing the gridded lidar data.  
  
  **Comments on parameters `t_win` and `x_win`**:  
  In some way, `t_win` and `x_win` are used to select the interpolation data points and indirectly define regions outside of which, data is not interpolated. So increasing those two parameters will fill in more gaps, but potentially increase the number of irrealistic points. This is slightly less true for `t_win`, but keep in mind that increasing `t_win` also slows down the interpolation process since more data is being used to create the interpolant.  

</details>


