function [ grid_data ] = fun_singlebeam_lidar_gridding( time , x_grid , raw_data , t_win , x_win )
% Function gridding (both time and space) the raw data from singlebeam lidar systems.
% This was initially written for processing the data collected at Duck with UAV-mounted lidar systems.
% It assumes minimal pre-processing for the input data, with the idea to make this function as generic as possible. 
% At minima, it needs a time (array of size Np x 1) and xyz (array of size Np x 3) to work properly.
%
% Inputs:
%   time     - time interpolation grid
%   x_grid   - x interpolation grid
%   raw_data - data structure containing de-noised raw data; minimum data fields: time (Np x 1); xyz (Np x 3)
%   t_win    - time window [s] within which data is used for the time interpolation 
%              e.g., if t_win = 0.15 s, then for each time ti of the grid, the algorithm will use data within 0.15 s 
%              both in past (ti-0.15) and future (ti+0.15). I suggest +- 1.5/original fs (i.e., 3 full scans are used)
%   x_win    - cross-shore space window [m] within which data needs to be present for keeping interpolated value (NaN is used otherwise)
%
% Outputs:
%   grid_data - a self-explanatory data structure containing the gridded lidar data
%
% Comments on parameters t_win and x_win:
%   In some way, t_win and x_win are used to select the interpolation data points and indirectly define regions outside of which, 
%   data is not interpolated. So increasing those two parameters will fill in more gaps, but potentially increase the number of irrealistic points.
%   This is slightly less true for t_win, but keep in mind that increasing t_win also slows down the interpolation process 
%   since more data is being used to create the interpolant.
%
% January 21, 2025
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Initialisation of gridded array
  grid_data.metadata  = ['Gridded lidar data starting on ',datestr(time(1))];
  grid_data.sf        = floor(1./(time(2)-time(1))/24/3600);
  grid_data.sf_info   = 'Sampling frequency of gridded data [Hz]';
  grid_data.time      = time;
  grid_data.time_info = 'UTC time [MATLAB format]';
  grid_data.dx        = x_grid(2)-x_grid(1);
  grid_data.dx_info   = 'Spatial resolution of gridded data [m]';
  grid_data.x         = x_grid;
  grid_data.x_info    = 'Cross-shore spatial grid [m]';
  grid_data.z         = nan(length(grid_data.time),length(grid_data.x));
  grid_data.z_info    = 'Interpolated elevation in the provided datum [m]';
  clear tbeg_vec tend_vec tbeg tend

  % Starting when we have the first data point and stopping at last
  first_step = min(raw_data.time); last_step = max(raw_data.time);
  itstart = find( grid_data.time > first_step , 1 , 'first' );
  itend   = find( grid_data.time < last_step , 1 , 'last' );

  % Making sure time indices are sorted; this will speed up algorithm
  if ~issorted(raw_data.time)
    [raw_data.time,isorted] = sort(raw_data.time);
    raw_data.xyz = raw_data.xyz(isorted,:);
  end

  % Main loop
  for tt = itstart:itend
    disp(['Number of timesteps left : ',num2str(length(grid_data.time)-tt)])
    % Looking for points within the correct time domain (typically a few time steps; imposed through "t_win")
    % Scarcity of points in the spatial domain is handled at a later stage, where interpolated points are removed if not enough data
%     iloc_t = find( abs(grid_data.time(tt) - raw_data.time ) <= t_win/24/3600 ); %original, much slower than by searching manually in sorted array, see below
    iloc_t = BinaryIntervalSearch(raw_data.time,grid_data.time(tt),t_win/24/3600)'; % taken from File-exchange

    % Interpolation, main part
    % Entering only when we have sufficient spatio-temporal data. Quality criteria hard-coded for now.
    if ~isempty(iloc_t) && length(iloc_t) > 5
      % Preparing interpolant
      interFunction = scatteredInterpolant(double(raw_data.time(iloc_t)),double(raw_data.xyz(iloc_t,1)),double(raw_data.xyz(iloc_t,3)),'linear','none');

      % Looping over spatial grid points
      % Interpolation only performed when sufficient data is present. Quality criteria hard-coded for now.
      for xx = 1:length(grid_data.x)
        % Looking for points within the correct time and spatial domains
        xloc = x_grid(xx);
        iloc_x = find( abs(raw_data.xyz(iloc_t,1)-xloc) < x_win );
        iloc   = iloc_t(iloc_x); %#ok<FNDSB> % Final indices

        % Interpolation only if linearisation is possible
        if and( ~isempty(iloc) , numel(iloc) > 2 )
          [Xq,Yq]            = ndgrid(time(tt),xloc); % Getting inputs in the correct dimensions
          grid_data.z(tt,xx) = interFunction(Xq,Yq);
        end
      end
    end
  end

  return
end
