% Workflow example for the gridding of multibeam lidar data. We start from a "pre-processed" (raw data, but organised such that it is read by the gridding function).
% The data originates from the BELS2022 experiments and corresponds to flight number #2 on the 12 September 2022 (1900-1930); lidar: Velodyne VLP32C, hovering position at x = 225 m.
% This script was prepared while working towards the paper:
% Seamless nearshore topo-bathymetry reconstruction from lidar scanners: a Proof-of-Concept based on a dedicated field experiment at Duck, NC
% by Martins Kévin, K.L. Brodie, J.W. Fiedler, A.M. O'Dea, N.J. Spore, R.L. Grenzeback, P.J. Dickhudt, S. Bak, O. de Viron and P. Bonneton 
% submitted to Coastal Engineering
%
% January 28, 2025
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all

% Matlab libraries
addpath(genpath('../'))

%% 1 - interpolating raw lidar data into time and space grids (full 4D approach)
% Interpolation grids
dx = 0.2;             % Spatial resolution
sf = 2;               % Sampling frequency of gridded timeframe
x_grid = [50:dx:300]; % Spatial resolution of gridded data

% Setting interpolation parameters
x_win = 1;    % Spatial window for using adjacents profiles
y_win = 1.5;  % Spatial window for using adjacents profiles
t_win = 0.16; % Temporal window in seconds. For t_win = 0.16 s and lidar rotations at 10 Hz, this means up to 3 scans can be used
              % Going higher is not recommended as it smooths waves face/front, while there are already enough data points for interpolation

% Time grid, fixed from the dune lidar schedule
time = [ datenum(2022,9,12,19,0,0) : 1/(sf*24*3600) : datenum(2022,9,12,19,30,0) ]';

% Dealing with existing files: do we want to re-update every pre-processed files?
update = false;

% Making sure we have the pre-processed data upon which the example is based
infilename = 'data/BELS_20220912_191123_flight_2_x=225.preprocessed.mat';
if and(~isfile(infilename),isunix) 
  disp('pre-processed file not found, downloading now...');
  system('pip install gdown && gdown https://drive.google.com/uc?id=13Qk01nyErT1LT-RtmG7dNQtoTl9gqpHI')
  system('mv 20220912_191123_flight_2_x=225.mat data/')
else
  error('pre-processed file not found, please download manually at: https://drive.google.com/uc?id=13Qk01nyErT1LT-RtmG7dNQtoTl9gqpHI');
end

% Gridding data
outfilename = ['data/BELS_20220912_1900_flight_2_x=225.dx=',num2str(dx),'m_',num2str(sf),'Hz_twin=',num2str(t_win),'s.mat'];
if ~isfile(outfilename) || update
  % Loading data
  raw_data = load( infilename );

  % Quick filtering of potentially abnomalous data
  idelete = find(or(or(raw_data.xyz(:,1) > 305,raw_data.xyz(:,1) < 170),or(raw_data.xyz(:,3) > 4,raw_data.xyz(:,3) < -2.5)));
  raw_data.beam_ID(idelete) = []; raw_data.time(idelete) = []; raw_data.xyz(idelete,:) = []; raw_data.intensity(idelete,:) = []; clear idelete

  % Gridding data
  y_grid    = x_grid*0 + 945;
  grid_data_1 = fun_multibeam_lidar_gridding( time , x_grid , y_grid , raw_data , t_win , x_win , y_win );

  % Saving gridded data
  save(outfilename,'-v7.3','-struct','grid_data_1');
else
  grid_data_1 = load(outfilename);
end
clear raw_data % freeing memory

% Checking dataset with basic information (here, we will work from x ~ 170 m to x ~ 250 m
fun_gridded_lidar_diagnostics( grid_data_1.sf , grid_data_1.x , grid_data_1.z , 1 );

%% 2 - interpolating raw lidar data into time and space grids (alongshore-uniform approach, faster but potentially less accurate)
% Here we use the singlebeam function, neglecting potential variations in y, which could be very advantageous for lab experiments.
% The raw data has to be restricted to a narrow swath of data, otherwise there is little speed-up gain.
% As we will see, it makes little differences in terms of final results here due to the normal propagation character of the waves that day.
% Interpolation grids
dx = 0.2;             % Spatial resolution
sf = 2;               % Sampling frequency of gridded timeframe
x_grid = [50:dx:300]; % Spatial resolution of gridded data

% Setting interpolation parameters
x_win = 1;    % Spatial window for using adjacents profiles
t_win = 0.16; % Temporal window in seconds. For t_win = 0.16 s and lidar rotations at 10 Hz, this means up to 3 scans can be used
              % Going higher is not recommended as it smooths waves face/front, while there are already enough data points for interpolation

% Time grid, fixed from the dune lidar schedule
time = [ datenum(2022,9,12,19,0,0) : 1/(sf*24*3600) : datenum(2022,9,12,19,30,0) ]';

% Dealing with existing files: do we want to re-update every pre-processed files?
update = false;

% Gridding data
infilename  = 'data/BELS_20220912_191123_flight_2_x=225.preprocessed.mat';
outfilename = ['data/BELS_20220912_1900_flight_2_x=225.dx=',num2str(dx),'m_',num2str(sf),'Hz_twin=',num2str(t_win),'s_noYdep.mat'];
if ~isfile(outfilename) || update
  % Loading data
  raw_data = load( infilename );

  % Quick filtering of potentially abnomalous data
  idelete = find(or(or(raw_data.xyz(:,1) > 305,raw_data.xyz(:,1) < 170),or(raw_data.xyz(:,3) > 4,raw_data.xyz(:,3) < -2.5)));
  idelete = cat(1,idelete,find(or(raw_data.xyz(:,2) > 946,raw_data.xyz(:,2) < 944))); % Narrowing it down
  raw_data.beam_ID(idelete) = []; raw_data.time(idelete) = []; raw_data.xyz(idelete,:) = []; raw_data.intensity(idelete,:) = []; clear idelete

  % Gridding data
  grid_data_2 = fun_singlebeam_lidar_gridding( time , x_grid , raw_data , t_win , x_win );

  % Saving gridded data
  save(outfilename,'-v7.3','-struct','grid_data_2');
else
  grid_data_2 = load(outfilename);
end
clear raw_data % freeing memory

% Checking dataset with basic information (here, we will work from x ~ 170 m to x ~ 250 m
fun_gridded_lidar_diagnostics( grid_data_2.sf , grid_data_2.x , grid_data_2.z , 1 );

%% Comparison
% Re loading data files in case you haven't executed the code above
dx = 0.2; sf = 2; t_win = 0.16; 
outfilename = ['data/BELS_20220912_1900_flight_2_x=225.dx=',num2str(dx),'m_',num2str(sf),'Hz_twin=',num2str(t_win),'s.mat'];
grid_data_1 = load(outfilename);
outfilename = ['data/BELS_20220912_1900_flight_2_x=225.dx=',num2str(dx),'m_',num2str(sf),'Hz_twin=',num2str(t_win),'s_noYdep.mat'];
grid_data_2 = load(outfilename);

% x = 220 m
[~,iloc] = nanmin(abs(grid_data_1.x-220));

% Timeseries
figure(1),
plot( grid_data_1.time, grid_data_1.z(:,iloc) ), hold on, box on, grid on
plot( grid_data_2.time, grid_data_2.z(:,iloc) )

% Spectral analysis
% Parameters for FFT and quality checks
nfft = 128*sf;  % Number of points used per FFT - i.e. block length
overlap = 75;   % Amount of overlap in % between blocks of data
wind = 'hann';  % Tapering: use hann or rectangular
thperNaN = 10;  % 10 to 15% is completely acceptable from first tests

% Re-organising data by block
[~,zeta_mat_1] = fun_prep_gappy_series_by_block( grid_data_1.time, grid_data_1.z(:,iloc) , nfft , overlap , thperNaN );
[~,zeta_mat_2] = fun_prep_gappy_series_by_block( grid_data_2.time, grid_data_2.z(:,iloc) , nfft , overlap , thperNaN );

% Computing PSD
psd_zeta_1 = fun_compute_spectrum_mat( zeta_mat_1 , sf , overlap , wind );
psd_zeta_2 = fun_compute_spectrum_mat( zeta_mat_2 , sf , overlap , wind );

% Spectra
scrsz = get(0,'ScreenSize'); fig1 = figure(2);
set(fig1,'Position',[550 600 scrsz(3)*0.15 scrsz(4)*0.25],'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 6 4],'color',[250,250,250]/255);
set(0,'defaultAxesFontSize',5), set(0,'defaultaxeslinewidth',0.5)

loglog( psd_zeta_1.f , psd_zeta_1.E ), hold on, box on, grid on
loglog( psd_zeta_2.f , psd_zeta_2.E )
% Axes
set(gca,'ytick',[10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1 100]); set(gca,'ylim',[10^-3 10])
set(gca,'xtick',[0.001 0.01 0.1 1 10]); set(gca,'xlim',[5*10^-3 3])
xlabel('f [Hz]'); ylabel('E(f) [m^2/Hz]');
% Add confidence levels
coef_ulimit = psd_zeta_1.CI(2); coef_llimit = psd_zeta_1.CI(1);
f_CI = 0.02; psd_CI = 0.025;
std_upper_limit = psd_CI*coef_ulimit;
std_lower_limit = psd_CI*coef_llimit;
errorbar( f_CI , psd_CI , abs(psd_CI-std_lower_limit) , abs(psd_CI-std_upper_limit) , 'ko' , 'markerfacecolor' , 'k' , 'markersize' , 2 , 'linewidth', 0.5); hold off
text( 0.013 , 0.08 , '95% C.I.','Fontsize',8); hold off
