# lidar-toolbox: a MATLAB set of libraries to handle nearshore wave data collected with lidars

This is part of **lidar-toolbox**: a MATLAB set of libraries to handle nearshore wave data collected with lidars. This sub-directory **gappy-series-preprocessing** gather useful functions for pre-processing gappy series of data so that spectral, cross-spectral and bispectral analyses can be applied to them (compatible with spectral-analysis library). This set of functions was initially created for working with ocean waves free surface elevation gappy signals obtained with lidars.  

<strong>Latest updates:</strong>  
<sub><sup>:arrow_forward:</sup></sub> *(Jan. 2025)*
gappy-series-preprocessing v1: first release of the library; example to be added separately, while reseasing entire lidBathy workflow.

<strong>Contact me:</strong>  
Kévin Martins  
CNRS researcher at UMR 7266 LIENSs, CNRS - La Rochelle University, France  
kevin.martins@cnrs.fr

---

## List of functions

<details>
  <summary>📄 fun_count_pNaNs.m</summary>  
  <br>  

  **Description**:  
  Function computing and returning the percentage of NaNs within a timeseries.

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `data`    | double | series of data (potentially) containing NaNs                     |

  **Outputs**:  
  &nbsp;&nbsp;Returns `pNaNs`, the percentage of NaNs within `data`.

</details>

<details>
  <summary>📄 fun_interp_series.m</summary>  
  <br>  

  **Description**:  
  Linear interpolation of a timeseries, with special treatment at the timeseries tails. More options could be later added, depending on the complexity required here, but for spectral analysis, sensitivity tests showed little influence so far.  

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `x`       | double | grid data (generally time, as mostly intended for timeseries) [whatever unit] |  
  | `data`    | double | series of data containing NaNs [whatever unit]; same length as `x`                                         |  
  | `tail_method` | int   | 1 - imposing first (last) non-NaN value at beginning (end) of series (preferred); or 2 - removing NaN values at beginning and end of series (NB: thus changing length of series) |

  **Outputs**:  
  &nbsp;&nbsp;`x`, the grid data, only useful when `tail_method = 2`  
  &nbsp;&nbsp;`data`, the interpolated series

</details>

<details>
  <summary>📄 fun_prep_gappy_series_by_block.m</summary>  
  <br>  

  **Description**:  
  Function restructuring a potentially gappy series into a matrix with overlapping blocks of length nfft. This is intended to prepare the timeseries for spectral analysis, using a threshold on the percentage of NaNs allowed per block, for quality control. This is particularly useful for lidar data, for instance, which can be naturally gappy, or to deal with large blocs of NaNs within a data series. We also keep track of time, if needed.  

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `time`    | double | grid data (generally time, as mostly intended for timeseries) [whatever unit] |  
  | `signal`  | int    | series of data potentially containing NaNs [whatever unit]; same length as `time`   |  
  | `nfft`    | int    | block length for the FFT (default = 256)                        |
  | `overlap` | int    | percentage overlap (typical is 50%, 75% optimises edof)        |
  | `thperNaN` | int   | maximal percentage of NaNs allowed within block of data |

  **Outputs**:  
  &nbsp;&nbsp; `time_mat`, time matrix corresponding to block of data.  
  &nbsp;&nbsp; `signal_mat`, data matrix of size ( nfft , number of blocks with less than thperNaN % of NaNs).  

</details>

<details>
  <summary>📄 fun_prep_gappy_series_by_block_xs.m</summary>  
  <br>  

  **Description**:  
  Function restructuring two potentially gappy series into matrices with overlapping blocks of length nfft. This is intended to prepare the timeseries for cross-spectral analysis, using a threshold on the percentage of NaNs allowed per block, for quality control. This is particularly useful for lidar data, for instance, which can be naturally gappy, or to deal with large blocs of NaNs within a data series. We also keep track of time, if needed.  

  **Inputs**:  

  | Name      | Type   | Description                                                      |
  |-----------|--------|------------------------------------------------------------------|
  | `time`    | double | grid data (generally time, as mostly intended for timeseries) [whatever unit] |  
  | `signal_1` | double | first series of data potentially containing NaNs [whatever unit]; same length as `time`   | 
  | `signal_2` | double | second series of data potentially containing NaNs [whatever unit]; same length as `time`  |  
  | `nfft`    | int    | block length for the FFT (default = 256)                        |
  | `overlap` | int    | percentage overlap (typical is 50%, 75% optimises edof)                     |
  | `thperNaN` | int   | maximal percentage of NaNs allowed within block of data |

  **Outputs**:  
  &nbsp;&nbsp; `time_mat`, time matrix corresponding to block of data.  
  &nbsp;&nbsp; `signal_1_mat`, data matrix from signal_1 of size ( nfft , number of blocks with less than thperNaN % of NaNs).  
  &nbsp;&nbsp; `signal_2_mat`, data matrix from signal_2 of size ( nfft , number of blocks with less than thperNaN % of NaNs).  

</details>

