function [ time_mat , signal_mat ] = fun_prep_gappy_series_by_block( time , signal , nfft , overlap , thperNaN )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function restructuring a potentially gappy series into a matrix with overlapping blocks of length nfft.
% This is intended to prepare the timeseries for spectral analysis, using a threshold on the percentage of NaNs allowed per block,
% for quality control. This is particularly useful for lidar data, for instance, which can be quite gappy, 
% or to deal with large blocs of NaNs within a data series. We also keep track of time, if needed.
%
% Inputs:
%   time       - time (unit not important, we send back what is input); size : length(time) x 1
%   signal     - series of (potentially gappy) data; size : length(time) x 1
%   nfft       - bloc length for the FFT [default = 256]
%   overlap    - percentage overlap (typical is 50%, 75% optimises edof)
%   thperNaN   - maximal percentage of NaNs allowed within block of data
%
% Outputs: 
%   time_mat   - time matrix corresponding to block of data
%   signal_mat - data matrix of size ( nfft , number of blocks with less than thperNaN % of NaNs)
%
% January 23, 2024
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --------------------- Various parameters -------------------------

  lt = size(signal,1); time = time(:);
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

  time_mat   = nan(nfft,nblock); 
  signal_mat = nan(nfft,nblock); 
  
  % Loop over blocks
  locseg = [1:nfft]'; % Initialising indices
  for kk = 1:nblock
    % Adjustment to locseg if this is the last block
    if (kk == nblock)
      locseg = [lt-nfft+1:lt]';
    end
    
    % Raw timeseries
    time_mat(:,kk) = time(locseg);
    signal_loc = signal(locseg);

    % Checking for NaNs
    if fun_count_pNaNs(signal_loc) < thperNaN
      % Dealing with NaNs (interpolation)
      [~,signal_loc] = fun_interp_series( [1:nfft]' , signal_loc , 1 );

      % Detrending
      signal_loc = detrend(signal_loc - nanmean(signal_loc));
    
      % Storing data in matrix
      signal_mat(:,kk) = signal_loc;

      % Marking block for keeping
      boolkeep = cat(1,boolkeep,kk);
    end
    
    % Indices for next block and updating count
    locseg = locseg + nadvance;
  end

  % We keep only marked blocks
  time_mat   = time_mat(:,boolkeep);
  signal_mat = signal_mat(:,boolkeep);
  
  return
end
