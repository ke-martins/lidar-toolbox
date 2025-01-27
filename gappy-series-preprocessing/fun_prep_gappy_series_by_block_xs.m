function [ time_mat , signal_1_mat , signal_2_mat ] = fun_prep_gappy_series_by_block_xs( time , signal_1 , signal_2 , nfft , overlap , thperNaN ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function restructuring two potentially gappy series into matrices with overlapping blocks of length nfft.
% This is intended to prepare the timeseries for cross-spectral analysis, using a threshold on the percentage of NaNs allowed per block,
% for quality control. This is particularly useful for lidar data, for instance, which can be quite gappy, 
% or to deal with large blocs of NaNs within a data series. We also keep track of time, if needed.
%
% Inputs:
%   time         - time (unit not important, we send back what is input); size : length(time) x 1
%   signal_1     - series 1 of (potentially gappy) data; size : length(time) x 1
%   signal_2     - series 2 of (potentially gappy) data; size : length(time) x 1
%   nfft         - bloc length for the FFT [default = 256]
%   overlap      - percentage overlap (typical is 50%, 75% optimises edof)
%   thperNaN     - maximal percentage of NaNs allowed within block of data
%
% Outputs: 
%   time_mat     - time matrix corresponding to block of data
%   signal_1_mat - data matrix from signal_1 of size ( nfft , number of blocks where both signal_1 and signal_2 have less than thperNaN % of NaNs)
%   signal_2_mat - data matrix from signal_2 of size ( nfft , number of blocks where both signal_1 and signal_2 have less than thperNaN % of NaNs)
%
% January 23, 2024
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --------------------- Various parameters -------------------------

  lt = size(signal_1,1); time = time(:);
  if (exist('nfft') ~= 1)         nfft = 256; end
  if (exist('overlap') ~= 1)   overlap = 50;  end
  if (exist('thperNaN') ~= 1) thperNaN = 10;  end

  nfft     = nfft - rem(nfft,2);
  overlap  = min(99,max(overlap,0));
  eadvance = fix(nfft * overlap / 100);
  nadvance = nfft - eadvance;
  nblock   = fix((lt - eadvance) / nadvance) + 1; % +1 for not throwing any data out
  boolkeep = [];
  

  % ---------------------- Initialization ------------------------

  time_mat     = nan(nfft,nblock); 
  signal_1_mat = nan(nfft,nblock); 
  signal_2_mat = nan(nfft,nblock); 
  
  % Loop over blocks
  locseg = [1:nfft]'; % Initialising indices
  for kk = 1:nblock
    % Adjustment to locseg if this is the last block
    if (kk == nblock)
      locseg = [lt-nfft+1:lt]';
    end
    
    % Raw timeseries
    time_mat(:,kk) = time(locseg);
    signal_1_loc   = signal_1(locseg);
    signal_2_loc   = signal_2(locseg);

    % Checking for NaNs
    if and(fun_count_pNaNs(signal_1_loc) < thperNaN,fun_count_pNaNs(signal_2_loc) < thperNaN)
      % Dealing with NaNs (interpolation)
      [~,signal_1_loc] = fun_interp_series( [1:nfft]' , signal_1_loc , 1 );
      [~,signal_2_loc] = fun_interp_series( [1:nfft]' , signal_2_loc , 1 );

      % Detrending
      signal_1_loc = detrend(signal_1_loc - nanmean(signal_1_loc));
      signal_2_loc = detrend(signal_2_loc - nanmean(signal_2_loc));
    
      % Storing data in matrix
      signal_1_mat(:,kk) = signal_1_loc;
      signal_2_mat(:,kk) = signal_2_loc;

      % Marking block for keeping
      boolkeep = cat(1,boolkeep,kk);
    end
    
    % Indices for next block and updating count
    locseg = locseg + nadvance;
  end

  % We keep only marked blocks
  time_mat     = time_mat(:,boolkeep);
  signal_1_mat = signal_1_mat(:,boolkeep);
  signal_2_mat = signal_2_mat(:,boolkeep);
end
