% Workflow example for the processing and preparation of lidar data prior to Boussinesq-based depth inversion
% This is flight number #1 on the 12 September 2022; lidar: Velodyne VLP32C, hovering position at x = 225 m
% This script comes as a supporting material for the paper:
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

%% 1 - interpolating raw lidar data into time and space grids
% Note: here we use a pre-processed data structure. This pre-processing is done elsewhere, as this is published work (e.g., see Fiedler et al., 2021)
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

% Gridding data
infilename  = '../examples/data/BELS_20220912_191123_flight_1_x=225.preprocessed.mat';
outfilename = ['../examples/data/BELS_20220912_1900_flight_1_x=225.dx=',num2str(dx),'m_',num2str(sf),'Hz_twin=',num2str(t_win),'s.mat'];
if ~isfile(outfilename) || update
  % Loading data
  raw_data = load( infilename );

  % Quick filtering of potentially abnomalous data
  idelete = find(or(or(raw_data.xyz(:,1) > 305,raw_data.xyz(:,1) < 170),or(raw_data.xyz(:,3) > 4,raw_data.xyz(:,3) < -2.5)));
  raw_data.beam_ID(idelete) = []; raw_data.time(idelete) = []; raw_data.xyz(idelete,:) = []; raw_data.intensity(idelete,:) = []; clear idelete

  % Gridding data
  y_grid    = x_grid*0 + 945;
  grid_data = fun_multibeam_lidar_gridding( time , x_grid , y_grid , raw_data , t_win , x_win , y_win );

  % Saving gridded data
  save(outfilename,'-v7.3','-struct','grid_data');
else
  grid_data = load(outfilename);
end
clear raw_data % freeing memory

% Checking dataset with basic information (here, we will work from x ~ 170 m to x ~ 250 m
fun_gridded_lidar_diagnostics( grid_data.sf , grid_data.x , grid_data.z , 1 );

% A bit of cleaning
clear infilename outfilename update x_win y_win t_win dx sf x_grid time

%% 2 - retrieving survey data along this transect
% Dealing with survey
survey = load('../examples/data/BELS_FRF_crawler+crab_survey_20220912_NAVD88.mat');

% Dealing with survey and water depth
interFunction = scatteredInterpolant([survey.x,survey.x]',[survey.y(1,:),survey.y(2,:)]',[survey.z(1,:),survey.z(2,:)]','linear','none');
grid_data.zb   = interFunction(grid_data.x,grid_data.y); clear interFunction
grid_data.zb_info = 'Interpolated seabed elevation [m] above NAVD88 datum'; 
clear survey
% figure(1), plot( grid_data.x , grid_data.zb ), hold on, plot(survey.x,survey.z(1,:),survey.x,survey.z(2,:))

%% 3 - initialisation of data structure for depth-inversion
% Parameters for analyses
p.thperNaN = 15;               % Maximal percentage of NaNs allowed within block of data
p.wind     = 'hann';           % Windowing applied to each block of data (only for cross-spectral analysis)
p.fLp_max  = 20;               % Maximal fraction of peak wavelength used to determine distance between sensors
p.fLp_min  = 8;                % Minimal fraction of peak wavelength used to determine distance between sensors
p.ixmean   = 0.2;              % Maximal distance between sensors used for getting surrounding information
p.overlap  = 75;               % Overlap in percentage
p.nfft     = 128*grid_data.sf; % For FFT
p.lt       = length(grid_data.time); eadvance = fix(p.nfft * p.overlap / 100);
p.nadvance = p.nfft - eadvance;
p.nblock   = fix((p.lt - eadvance) / p.nadvance) + 1; % +1 for not throwing any data out
clear eadvance

% Preparing depth-inversion array
data.p               = p;
data.p_info          = 'List of parameters for performing the analysis';
data.x               = [210,225,240]; % Performing analysis only on a subset, since we use data from SIO lidar only
data.x_info          = grid_data.x_info; 
data.mwl             = interp1(grid_data.x,nanmean(grid_data.z),data.x);
data.mwl_info        = 'Mean water level [m above NAVD88]';
data.hm              = interp1(grid_data.x,nanmean(grid_data.z)-grid_data.zb,data.x);
data.hm_info         = grid_data.x_info;
data.c_xcor          = nan(1,length(data.x));
data.c_xcor_info     = 'Bulk celerity [m/s] from simple cross-correlation (mean of several combination)';
data.c_xcor_std      = nan(1,length(data.x));
data.c_xcor_std_info = 'Std of bulk celerity [m/s] from simple cross-correlations';
data.f               = [0:p.nfft/2]'/p.nfft*grid_data.sf;
data.f_info          = 'Frequencies [Hz]';
data.E               = nan((p.nfft)/2 + 1,length(data.x));
data.E_info          = 'Energy density spectrum [m^2/Hz]';
data.E_CI            = nan(2,length(data.x));
data.E_CI_info       = 'Confidence interval';
% Bispectrum part
data.bis.f           = [-p.nfft/2:p.nfft/2]'/p.nfft*grid_data.sf;
data.bis.f_info      = 'Frequencies [Hz]';
data.bis.P           = nan(p.nfft+1,length(data.x));
data.bis.P_info      = 'Power spectrum [m^2] - no tapering';
data.bis.B           = nan(p.nfft+1,p.nfft+1,length(data.x));
data.bis.B_info      = 'Power bispectrum [m^3] - no tapering';
data.bis.B_std       = nan(p.nfft+1,p.nfft+1,length(data.x));
data.bis.B_std_info  = 'Power bispectrum standard deviation [m^3]';
data.bis.k_rms       = nan(p.nfft+1,length(data.x));
data.bis.k_rms_info  = 'Predicted wavenumbers [m^-1] using Boussinesq theory';
% Cross-spectrum part
data.xspe.f          = [0:p.nfft/2]'/p.nfft*grid_data.sf;
data.xspe.f_info     = 'Frequencies [Hz]';
data.xspe.k_rms      = nan((p.nfft)/2 + 1,length(data.x));
data.xspe.k_rms_info = 'Mean wavenumber [m^-1] from cross-spectral analysis';
data.xspe.k_rms_std  = nan((p.nfft)/2 + 1,length(data.x));
data.xspe.k_rms_std_info = 'Wavenumber std [m^-1] from cross-spectral analysis';
data.xspe.k_rms_CI   = nan((p.nfft)/2 + 1,length(data.x));
data.xspe.k_rms_CI_info = 'Wavenumber confindence interval [m^-1] from cross-spectral analysis';
data.xspe.coh        = nan((p.nfft)/2 + 1,length(data.x));
data.xspe.coh_info   = 'Coherence [-] from cross-spectral analysis';

%% 4 - performing bispectral analysis
% Loop over cross-shore positions
for ix = 1:numel(data.x)
  disp(['Bispectral analysis at x = ',num2str(data.x(ix)),' m. ',num2str(numel(data.x)-ix),' positions left.'])

  % Looking where to work in original dataset
  [~,xx] = nanmin(abs(grid_data.x-data.x(ix)));% grid_data.x(xx)

  % Adjustments around sensors, used for increasing stability
  ildx = round(p.ixmean/grid_data.dx);
  iidx = -ildx:ildx; clear ildx

  % Gathering blocks of data around location into a matrix
  % Data blocks neet to meet a certain criteria based on NaN content
  zeta_mat = [];
  for ii = iidx
    [~,zeta] = fun_prep_gappy_series_by_block( grid_data.time , grid_data.z(:,xx+ii) , p.nfft , p.overlap , p.thperNaN );
    zeta_mat = cat(2,zeta_mat,zeta);
  end

  % Continuing only if we have enough data
  if size(zeta_mat,2) > 4
    % Computing spectrum and bispectrum
    bis_rect = fun_compute_bispectrum_mat( zeta_mat , grid_data.sf , p.overlap , 'rectangular' );
    psd_zeta = fun_compute_spectrum_mat( zeta_mat , grid_data.sf , p.overlap , p.wind );
    
    % Storing bispectrum products
    data.bis.B(:,:,ix)     = bis_rect.B;
    data.bis.B_std(:,:,ix) = bis_rect.B_std;
    data.bis.P(:,ix)       = bis_rect.P;
    data.E(:,ix)           = psd_zeta.E;
    data.E_CI(:,ix)        = psd_zeta.CI;
  end
end
clear ix ii ildx iidx xx zeta zeta_mat bis_rect psd_zeta

%% 5 - performing cross-spectral analysis
% Estimating the peak wave frequency
[~,xx] = nanmin(abs(grid_data.x-240));

% Adjustments around sensors, used for increasing stability
ildx = round(p.ixmean/grid_data.dx);
iidx = -ildx:ildx;

% Gathering blocks of data around location
zeta_mat = [];
for ii = iidx
  [~,zeta] = fun_prep_gappy_series_by_block( grid_data.time , grid_data.z(:,xx+ii) , p.nfft , p.overlap , p.thperNaN );
  zeta_mat = cat(2,zeta_mat,zeta);
end

% Computing PSD
psd_zeta = fun_compute_spectrum_mat( zeta_mat , grid_data.sf , p.overlap , p.wind );

% Plot
figure(2), loglog( psd_zeta.f , psd_zeta.E , 'k' )
fp = 0.0625;
clear ii xx iidx ildx zeta zeta_mat


%%%%%%%%%%%%%%%%%%%%%%
% Preparation of data
% Here we cheat a little bit by computing an approximate celerity, but we could measure it by cross-correlation...
hm  = max(nanmean(grid_data.z)-grid_data.zb,0);
Lfp = sqrt(9.81*hm)/fp; clear hm

% Looking for the maximum of possible combinations based on estimated wavelength (just for initialisation)
imax   = fix(min(max(1,p.fLp_max*Lfp/2/100),15)/grid_data.dx);
imin   = fix(max(min(10,p.fLp_min*Lfp/2/100),2)/grid_data.dx);
ncombi = max(imax)-min(imin);

% Max number of data blocks, useful to initialize some arrays
ildx = round(p.ixmean/grid_data.dx);
iidx = -ildx:ildx; clear ildx
n  = fix((100-p.overlap)/100*p.nfft);
Nb = numel(iidx)*floor((numel(grid_data.time)-p.nfft)/n+1); clear n


%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the cross-spectral arrays and computation
% Here we will store all data, from all possible combination of gauges around the location where we want to estimate phase velocity spectra
cpsd_data.f          =  [0:p.nfft/2]'/p.nfft*grid_data.sf;
cpsd_data.f_info     = 'Frequency [Hz] (one-sided)';
cpsd_data.df         = abs(cpsd_data.f(2)-cpsd_data.f(1));
cpsd_data.df_info    = 'Frequency resolution [Hz]';
cpsd_data.x          = data.x;
cpsd_data.x_info     = 'Cross-shore position [m]';
cpsd_data.dx         = nan(ncombi,numel(data.x));
cpsd_data.dx_info    = 'Distance between gauges [m] for the specific combination';
cpsd_data.cpsd       = nan(numel(cpsd_data.f),ncombi,numel(data.x));
cpsd_data.cpsd_info  = 'Cross-spectrum between two signals [-]';
cpsd_data.phase      = nan(numel(cpsd_data.f),ncombi,numel(data.x));
cpsd_data.phase_info = 'Phase spectrum between two signals [-]';
cpsd_data.coh        = nan(numel(cpsd_data.f),ncombi,numel(data.x));
cpsd_data.coh_info   = 'Coherence spectrum between two signals [-]';
cpsd_data.k_rms      = nan(numel(cpsd_data.f),ncombi,numel(data.x));
cpsd_data.k_rms_info = 'Root-mean square wavenumber [m^-1]';
cpsd_data.c_bulk     = nan(ncombi,numel(data.x));
cpsd_data.c_bulk_info = 'Mean wave propagation speed [m/s]';
cpsd_data.c_bulk_std = nan(ncombi,numel(data.x));
cpsd_data.c_bulk_std_info = 'Standard deviation associated with mean wave propagation speeds [m/s]';


%%%%%%%%%%%%%%%%%%%%%%
% Main code
% We work in a centered manner: all correlations are performed around a main location
for ix = 1:numel(data.x)
  disp(['Cross-spectral analysis at x = ',num2str(data.x(ix)),' m. ',num2str(numel(data.x)-ix),' positions left.'])

  % Looking where to work in original dataset
  [~,xx] = nanmin(abs(grid_data.x-data.x(ix)));

  % Quick filter on local water depth
  if nanmean(grid_data.z(:,xx))-grid_data.zb(xx) < 0.4
    disp(['Not deep enough, skipping position.'])
  end

  % Searching for all potential combinations
  combilist = imin(xx):imax(xx); 
  combi = [-combilist',combilist'];

  % Computing distance between gauges
  dx_loc = nan(1,numel(combilist));
  for cc = 1:numel(combilist)
    dx_loc(cc) = abs(grid_data.x(xx+combi(cc,1))-grid_data.x(xx+combi(cc,2)));
  end
  cpsd_data.dx(1:numel(dx_loc),ix) = dx_loc;

  % Computing cross-spectral products and mean propagation speed
  c_time   = nan(numel(combilist),Nb);
  c_bulk   = nan(numel(combilist),Nb);
  Nb_comb  = 0; % Effective number of combination used (i.e. with enough data)
  Nb_block = nan(1,numel(combilist)); % Number of blocks per combination used (needed to compute k_rms_CI)
  for cc = 1:numel(combilist)
    if and( xx+combi(cc,1) >= 1 , xx+combi(cc,2) <= numel(grid_data.x) )
      % Gathering blocks of data around location
      time_mat = []; zeta_l_mat = []; zeta_s_mat = [];
      for ii = iidx
        [time_loc,zeta_l,zeta_s] = fun_prep_gappy_series_by_block_xs( grid_data.time , grid_data.z(:,xx+combi(cc,1)+ii) , grid_data.z(:,xx+combi(cc,2)+ii) , ...
                                                                      p.nfft , p.overlap , p.thperNaN );
        time_mat   = cat(2,time_mat,time_loc);
        zeta_s_mat = cat(2,zeta_s_mat,zeta_s);
        zeta_l_mat = cat(2,zeta_l_mat,zeta_l);
      end

      % Continuing only if we have enough data
      if size(zeta_s_mat,2) > 4
        % Counting combination
        Nb_comb = Nb_comb + 1; Nb_block(cc) = size(zeta_s_mat,2);

        % Computing cross-spectrum, its coherence and phase
        cspec_loc = fun_compute_cross_spectrum_mat( zeta_l_mat , zeta_s_mat , grid_data.sf , p.wind );

        % Storing local stuff
        cpsd_data.cpsd(:,cc,ix)  = cspec_loc.cpsd;
        cpsd_data.coh(:,cc,ix)   = ((abs(cspec_loc.cpsd).^2)./(cspec_loc.Exx.*cspec_loc.Eyy)).^0.5;
        cpsd_data.phase(:,cc,ix) = atan2(-imag(cspec_loc.cpsd),real(cspec_loc.cpsd));
        cpsd_data.k_rms(:,cc,ix) = unwrap(cpsd_data.phase(:,cc,ix))/dx_loc(cc);

        % Computing mean bulk celerity
        for bb = 1:size(zeta_s_mat,2)
          c_time(cc,bb) = time_mat(bb);
          c_bulk(cc,bb) = abs(fun_compute_c_from_xcorr( zeta_l_mat(:,bb)' , zeta_s_mat(:,bb)' , grid_data.sf , dx_loc(cc) , 7.5 ));
        end

        % Storing bulk celerity data
        cpsd_data.c_bulk(cc,ix)     = nanmean(c_bulk(cc,:));
        cpsd_data.c_bulk_std(cc,ix) = nanstd(c_bulk(cc,:));
      end
    end
  end

  % Working out and saving data only if there were enough data
  if Nb_comb > 1
    % Bulk estimates
    data.c_xcor(ix)     = nanmean(cpsd_data.c_bulk(:,ix));
    data.c_xcor_std(ix) = nanmean(cpsd_data.c_bulk_std(:,ix));

    % First filter based on bulk celerity
    c_rms = 2*pi*cpsd_data.f ./ squeeze(cpsd_data.k_rms(:,:,ix));
    [idel,jdel] = find( or( c_rms > 2*data.c_xcor(ix) , c_rms < 0 ) );
    cpsd_data.k_rms(idel,jdel,ix)  = NaN;

    % Filtering outliers based on weird wave phase speeds
    k_rms = squeeze(cpsd_data.k_rms(:,:,ix)); 
    c_rms = 2*pi*cpsd_data.f ./ k_rms; %figure, plot( cpsd_data.f , c_rms)
    crms_mean = nanmean(c_rms); crms_std = nanstd(c_rms);
    [idel,jdel] = find( abs(c_rms-crms_mean) > 2*crms_std );
    cpsd_data.k_rms(idel,jdel,ix)  = NaN;

    % Working out variability and storing results
    % Spectral estimates
    data.xspe.k_rms(:,ix)     = nanmean(cpsd_data.k_rms(:,:,ix),2);
    data.xspe.k_rms_std(:,ix) = nanstd(cpsd_data.k_rms(:,:,ix),[],2); % this is useful to check that our theoretical CI is actually much larger than the variability of the cross-spectral estimate...
    data.xspe.coh(:,ix)       = nanmean(cpsd_data.coh(:,:,ix),2);
    var_k_rms = 1/(2*nanmean(Nb_block)) * nanmean((squeeze(cpsd_data.coh(:,1:numel(dx_loc),ix)).^(-2) - 1) ./ repmat(dx_loc,numel(data.xspe.f),1),2); %figure, plot( data.xspe.f , var_k_rms , 'k*' )
    data.xspe.k_rms_CI(:,ix)  = 2*sqrt(var_k_rms);
  end

  clear ii iii cc xx bb idel jdel c_rms crms_mean crms_std k_rms zeta_s_mat zeta_l_mat zeta_s zeta_l cspec_loc dx_loc time_mat time_loc var_k_rms Nb_comb Nb_block
end

%% 6 - depth-inversion and Monte Carlo approach
% Size of population for Monte-Carlo procedure
Npop = 5000;

% Initialising bathy data structure
bathy_data.x           = data.x;
bathy_data.x_info      = 'Cross-shore position [m]';
bathy_data.fp          = [fp,fp,fp];
bathy_data.fp_info     = 'Peak frequency [Hz] during flight';
bathy_data.mwl         = nan(1,numel(data.x));
bathy_data.mwl_info    = 'Mean water elevation [m above NAVD88]';
bathy_data.hm          = nan(1,numel(data.x));
bathy_data.hm_info     = 'Depth estimates [m] from lidar data and CRAB surveys';
bathy_data.h_krms      = nan(Npop,numel(data.x));
bathy_data.h_krms_info = 'Depth estimates [m] from fully-spectral Boussinesq approximation - using krms';
bathy_data.h_L         = nan(Npop,numel(data.x));
bathy_data.h_L_info    = 'Depth estimates [m] from fully-spectral linear theory';

% Looping over cross-shore locations
for xx = 1:numel(data.x)
  % MWL
  bathy_data.mwl(xx) = data.mwl(xx);

  % Mean water depth deduced from the two latter variables
  zb = interp1(grid_data.x, grid_data.zb, data.x(xx));
  bathy_data.hm(xx) = bathy_data.mwl(xx) - zb; clear zb

  % Parameters
  fp     = bathy_data.fp(xx);
  fmin   = 0.8*fp;
  fmax   = 0.25; % Frequency limits for minimisation
  ifreqs = find( and(data.bis.f>fmin,data.bis.f<fmax) );

  % Main loop of the Monte-Carlo approach to estimate depth uncertainty
  for nn = 1:Npop
    % Bispectrum, modified randomly within their theoretical 95% C.I. (uniform distribution)
    B_loc = squeeze(data.bis.B(:,:,xx));
    B_loc = squeeze(data.bis.B(:,:,xx)) + (- 1 + 2*rand(size(B_loc))).*squeeze(data.bis.B_std(:,:,xx));

    % Compute bispectral terms gamma_f and gamma_a (see Martins et al., 2023, for more details)
    [gamma_f_loc,~,gamma_a_loc] = fun_compute_krms_terms( data.bis.f , data.bis.P(:,xx) , B_loc );
    freqs   = data.bis.f(ifreqs);
    gamma_f = gamma_f_loc(ifreqs);
    gamma_a = gamma_a_loc(ifreqs);

    % Getting cross-spectrum observations on bispectrum frequency array (here similar df, but shifted indices)
    coh_obs = nan(size(freqs));
    k_rms   = nan(size(freqs));
    dk_rms  = nan(size(freqs));
    for ff = 1:numel(freqs)
      [~,iff] = nanmin(abs(freqs(ff)-data.xspe.f));
      coh_obs(ff) = data.xspe.coh(iff,xx);
      k_rms(ff)   = data.xspe.k_rms(iff,xx);
      dk_rms(ff)  = data.xspe.k_rms_CI(iff,xx);
    end
    clear ff iff

    % Reducing to non-NaN
    ids     = find(~isnan(k_rms));
    freqs   = freqs(ids);
    k_rms   = k_rms(ids);
    dk_rms  = dk_rms(ids);
    coh_obs = coh_obs(ids);
    gamma_f = gamma_f(ids);
    gamma_a = gamma_a(ids);
    clear ids

    % Observed wavenumbers, modified randomly within their theoretical 95% C.I. (uniform distribution)
    kappa_loc = k_rms + (- 1 + 2*rand(length(k_rms),1)).*dk_rms;

    % Boussinesq approach
    kappa_rms_fun = @(w,f,a,hh,gg) w./sqrt(gg*hh).*sqrt( 1 + hh.*f - 1./hh*a); % target function
    minimized_fun = @(h) sum(coh_obs.*abs(kappa_loc - kappa_rms_fun(2*pi*freqs,gamma_f,gamma_a,h,9.81)).^1); % norm on L1
    bathy_data.h_krms(nn,xx) = fminbnd(minimized_fun,0.3,6);

    % Linear estimate
    minimized_fun = @(h) sum(coh_obs.*abs(h - atanh((2*pi*freqs).^2./kappa_loc/9.81)./kappa_loc).^1); % norm on L1
    bathy_data.h_L(nn,xx) = fminbnd(minimized_fun,0.3,6);
  end
end

% Making sure we have our 95% C.I. within +-2std from the median
for xx = 1:numel(data.x)
  iset = find( abs( bathy_data.h_krms(:,xx) - nanmedian(bathy_data.h_krms(:,xx) ) ) < 2*nanstd(bathy_data.h_krms(:,xx) ) ); abs(95-numel(iset)/Npop*100) %error in %
  iset = find( abs( bathy_data.h_L(:,xx)  - nanmedian(bathy_data.h_L(:,xx) ) ) < 2*nanstd(bathy_data.h_L(:,xx) ) ); abs(95-numel(iset)/Npop*100) 
end

%% 7 - plotting results
scrsz = get(0,'ScreenSize'); fig1 = figure(3);
set(fig1,'Position',[550 600 scrsz(3)*0.25 scrsz(4)*0.5],'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 25 16]);
set(0,'defaultAxesFontSize',9), set(0,'defaultaxeslinewidth',0.5)

% PSD - x = 210 m
h(1) = subplot(3,3,1);
semilogy( 0  , NaN , 'color' , [0,0.5,1] , 'linewidth' , 1. ); hold on; grid on; box on;
hmp = fill( [ fmin fmax fmax fmin ] , [ 10^-6 10^-6 100 100 ] , [0.9,0.9,0.9] ); set(hmp,'EdgeColor',[0.9,0.9,0.9]);
loglog( data.f, data.E(:,1) , 'k' , 'linewidth' , 1. );
% xlabel('$f$ [Hz]','Interpreter','latex','Fontsize',11); 
ylabel('$E(f)$ [m$^2$/Hz]','Interpreter','latex','Fontsize',12); 
text(0.015,3,'a)','Fontsize',10); text(0.085,3,'x = 210 m','Fontsize',10);

% PSD - x = 225 m
h(2) = subplot(3,3,2);
semilogy( 0  , NaN , 'color' , [0,0.5,1] , 'linewidth' , 1. ); hold on; grid on; box on;
hmp = fill( [ fmin fmax fmax fmin ] , [ 10^-6 10^-6 100 100 ] , [0.9,0.9,0.9] ); set(hmp,'EdgeColor',[0.9,0.9,0.9]);
loglog( data.f, data.E(:,2) , 'k' , 'linewidth' , 1. );
text(0.015,3,'b)','Fontsize',10); text(0.085,3,'x = 225 m','Fontsize',10);

% PSD - x = 240 m
h(3) = subplot(3,3,3);
semilogy( 0  , NaN , 'color' , [0,0.5,1] , 'linewidth' , 1. ); hold on; grid on; box on;
hmp = fill( [ fmin fmax fmax fmin ] , [ 10^-6 10^-6 100 100 ] , [0.9,0.9,0.9] ); set(hmp,'EdgeColor',[0.9,0.9,0.9]);
loglog( data.f, data.E(:,3) , 'k' , 'linewidth' , 1. ); 
text(0.015,3,'c)','Fontsize',10); text(0.085,3,'x = 240 m','Fontsize',10);

% Axes
for hh = 1:3
  subplot(h(hh))
  set(gca,'ytick',[10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1 100]); set(gca,'ylim',[10^-4 10])
  set(gca,'Xtick',[0 0.05 0.1 0.2:0.1:1],'Xlim',[0 0.4])
  set(gca, 'layer', 'top'); set(gca,'TickDir','out','FontSize',9);
  % Confidence levels
  coef_ulimit = data.E_CI(2,hh); coef_llimit = data.E_CI(1,hh);
  f_CI = 0.15; psd_CI = 0.004;
  std_upper_limit = psd_CI*coef_ulimit;
  std_lower_limit = psd_CI*coef_llimit;
  errorbar( f_CI , psd_CI , abs(psd_CI-std_lower_limit) , abs(psd_CI-std_upper_limit) , 'ko' , 'markerfacecolor' , 'k' , 'markersize' , 3 , 'linewidth', 0.5); hold off
  txt_lb = text( 0.125 , 0.0105 , '95% C.I.','Fontsize',8); hold off; set(txt_lb, 'layer', 'front')
end

h(4) = subplot(3,3,4); pl = nan(1,4);
hmp = fill( [ fmin fmax fmax fmin ] , [ 0 0 8 8 ] , [0.9,0.9,0.9] ); hold on, grid on, box on; set(hmp,'EdgeColor',[0.9,0.9,0.9]);
% Explicit approximation for k_L based on Guo
k_L = (2*pi*data.f).^2/9.81 .* (1-exp(-((2*pi*data.f)*sqrt(data.hm(1)/9.81)).^(5/2))).^(-2/5);
c_L = 2*pi*data.f ./ k_L; c_L(1) = c_L(2);
pl(1) = loglog( data.f, c_L , 'r' , 'linewidth' , 1 );
pl(4) = plot( data.xspe.f, 2*pi*data.xspe.f ./ data.xspe.k_rms(:,1) ,'ko' , 'markersize' , 4.);
pl(3) = plot( data.bis.f, 2*pi*data.bis.f ./ fun_compute_krms( median(bathy_data.h_krms(:,1)) , data.bis.f , data.bis.P(:,1) , data.bis.B(:,:,1) ) , '*' , 'color' , [0.8,0.5,0.1], 'markersize' , 4.5 , 'linewidth' , 0.4 );
pl(2) = plot( data.bis.f, 2*pi*data.bis.f ./ fun_compute_krms( data.hm(1) , data.bis.f , data.bis.P(:,1) , data.bis.B(:,:,1) ) , 'k*', 'markersize' , 4.5 , 'linewidth' , 0.4 );
set(gca,'Xtick',[0 0.05 0.1 0.2:0.1:1],'Xlim',[0 0.5],'Ylim',[2 6])
% Legend
ylabel('$c$ [m/s]','Interpreter','latex','Fontsize',12);
leg = legend(pl,'$c_{L}$','$c_{rms}$ - survey','$c_{rms}$ - fit','$c_{obs}$','Location','South','Interpreter','Latex','Orientation','Horizontal','Fontsize',10,'NumColumns',2); legend boxoff
leg.ItemTokenSize = [14,14,14]; text(0.02,6.5,'d)','Fontsize',10);

h(5) = subplot(3,3,5); pl = nan(1,4);
hmp = fill( [ fmin fmax fmax fmin ] , [ 0 0 8 8 ] , [0.9,0.9,0.9] ); hold on, grid on, box on; set(hmp,'EdgeColor',[0.9,0.9,0.9]);
% Explicit approximation for k_L based on Guo
k_L = (2*pi*data.f).^2/9.81 .* (1-exp(-((2*pi*data.f)*sqrt(data.hm(2)/9.81)).^(5/2))).^(-2/5);
c_L = 2*pi*data.f ./ k_L; c_L(1) = c_L(2);
pl(1) = loglog( data.f, c_L , 'r' , 'linewidth' , 1 );
pl(4) = plot( data.xspe.f, 2*pi*data.xspe.f ./ data.xspe.k_rms(:,2) ,'ko' , 'markersize' , 4.);
pl(3) = plot( data.bis.f, 2*pi*data.bis.f ./ fun_compute_krms( median(bathy_data.h_krms(:,2)) , data.bis.f , data.bis.P(:,2) , data.bis.B(:,:,2) ) , '*' , 'color' , [0.8,0.5,0.1], 'markersize' , 4.5 , 'linewidth' , 0.4 );
pl(2) = plot( data.bis.f, 2*pi*data.bis.f ./ fun_compute_krms( data.hm(2) , data.bis.f , data.bis.P(:,2) , data.bis.B(:,:,2) ) , 'k*', 'markersize' , 4.5 , 'linewidth' , 0.4 );
set(gca,'Xtick',[0 0.05 0.1 0.2:0.1:1],'Xlim',[0 0.5],'Ylim',[2 6])
% Legend
leg = legend(pl,'$c_{L}$','$c_{rms}$ - survey','$c_{rms}$ - fit','$c_{obs}$','Location','South','Interpreter','Latex','Orientation','Horizontal','Fontsize',10,'NumColumns',2); legend boxoff
leg.ItemTokenSize = [14,14,14]; text(0.02,6.5,'e)','Fontsize',10);

h(6) = subplot(3,3,6); pl = nan(1,4);
hmp = fill( [ fmin fmax fmax fmin ] , [ 0 0 8 8 ] , [0.9,0.9,0.9] ); hold on, grid on, box on; set(hmp,'EdgeColor',[0.9,0.9,0.9]);
% Explicit approximation for k_L based on Guo
k_L = (2*pi*data.f).^2/9.81 .* (1-exp(-((2*pi*data.f)*sqrt(data.hm(3)/9.81)).^(5/2))).^(-2/5);
c_L = 2*pi*data.f ./ k_L; c_L(1) = c_L(2);
pl(1) = loglog( data.f, c_L , 'r' , 'linewidth' , 1 );
pl(4) = plot( data.xspe.f, 2*pi*data.xspe.f ./ data.xspe.k_rms(:,3) ,'ko' , 'markersize' , 4.);
pl(3) = plot( data.bis.f, 2*pi*data.bis.f ./ fun_compute_krms( median(bathy_data.h_krms(:,3)) , data.bis.f , data.bis.P(:,3) , data.bis.B(:,:,3) ) , '*' , 'color' , [0.8,0.5,0.1], 'markersize' , 4.5 , 'linewidth' , 0.4 );
pl(2) = plot( data.bis.f, 2*pi*data.bis.f ./ fun_compute_krms( data.hm(3) , data.bis.f , data.bis.P(:,3) , data.bis.B(:,:,3) ) , 'k*', 'markersize' , 4.5 , 'linewidth' , 0.4 );
set(gca,'Xtick',[0 0.05 0.1 0.2:0.1:1],'Xlim',[0 0.5],'Ylim',[2 6])
% Legend
leg = legend(pl,'$c_{L}$','$c_{rms}$ - survey','$c_{rms}$ - fit','$c_{obs}$','Location','South','Interpreter','Latex','Orientation','Horizontal','Fontsize',10,'NumColumns',2); legend boxoff
leg.ItemTokenSize = [14,14,14]; text(0.02,6.5,'f)','Fontsize',10);

% Axes
for hh = 4:6
  subplot(h(hh))
  xlabel('$f$ [Hz]','Interpreter','latex','Fontsize',11); set(gca, 'layer', 'top');
  set(gca,'Xtick',[0 0.05 0.1 0.2:0.1:1],'Xlim',[0 0.4],'Ylim',[2 7])
  set(gca, 'layer', 'top'); set(gca,'TickDir','out','FontSize',9);
end

h(7) = subplot(3,3,7);
df = 0.005; bin_axis = [0:df:4];
pl_L   = hist(bathy_data.h_L(:,1)/data.hm(1),bin_axis);
pl_rms = hist(bathy_data.h_krms(:,1)/data.hm(1),bin_axis);
ibeg = find(pl_L>0,1,'first'); iend = find(pl_L>0,1,'last');
plot( bin_axis(ibeg-1:iend+1)+df/2 , pl_L(ibeg-1:iend+1)/sum(pl_L(ibeg-1:iend+1)) , 'r' , 'linewidth' , 1 ); hold on, grid on, box on;
ibeg = find(pl_rms>0,1,'first'); iend = find(pl_rms>0,1,'last'); 
plot( bin_axis(ibeg-1:iend+1)+df/2 , pl_rms(ibeg-1:iend+1)/sum(pl_rms(ibeg-1:iend+1)) , 'k' , 'linewidth' , 1 )
% set(pl_rms,'edgecolor','k','linewidth',1)
% plot(data.hm,0.15,'s','markersize',4,'markerfacecolor',[0.8,0.5,0.1],'linewidth',1)
plot(median(bathy_data.h_krms(:,1)/data.hm(1)),0.015,'ko','markersize',4,'markerfacecolor','k')
pl_rms = errorbar(median(bathy_data.h_krms(:,1)/data.hm(1)),0.015,2*std(bathy_data.h_krms(:,1)/data.hm(1)),'horizontal','k','linewidth',0.75);
text(median(bathy_data.h_krms(:,1)/data.hm(1))-0.05,0.0098,'median','color','k','FontSize',8)
text(median(bathy_data.h_krms(:,1)/data.hm(1))-0.035,0.004,'\pm 2\sigma_h','color','k','FontSize',8)
plot(median(bathy_data.h_L(:,1)/data.hm(1)),0.015,'ro','markersize',4,'markerfacecolor','r')
pl_rms = errorbar(median(bathy_data.h_L(:,1)/data.hm(1)),0.015,2*std(bathy_data.h_L(:,1)/data.hm(1)),'horizontal','r','linewidth',0.75);
text(median(bathy_data.h_L(:,1)/data.hm(1))-0.06,0.0098,'median','color','r','FontSize',8)
text(median(bathy_data.h_L(:,1)/data.hm(1))-0.045,0.004,'\pm 2\sigma_h','color','r','FontSize',8)
set(gca,'Xtick',[0:0.1:4],'Xlim',[0.85 1.85],'Ylim',[0 0.06],'Ytick',[0:0.01:0.4]); set(gca, 'layer', 'top'), text(0.895,0.0525,'g)','Fontsize',10); hold off
% y Axis
ylabel('$p$ [-]','Interpreter','latex','Fontsize',12);
% Statistics
text(1.02,0.0535,['Surveyed {\it h} = ',num2str(fix(data.hm(1)*100)/100),' m'],'color',[0.6,0.6,0.6],'FontSize',9)
text(0.887,0.042,['Median = ',num2str(fix(median(bathy_data.h_krms(:,1))*100)/100),' m'],'FontSize',8)
text(0.89,0.041-0.0065,['\sigma_h = ',num2str(fix(std(bathy_data.h_krms(:,1))*100)/100),' m'],'FontSize',8)
text(1.417,0.044,['Median = ',num2str(fix(median(bathy_data.h_L(:,1))*100)/100),' m'],'Color','r','FontSize',8)
text(1.42,0.043-0.0065,['\sigma_h = ',num2str(fix(std(bathy_data.h_L(:,1))*100)/100),' m'],'Color','r','FontSize',8)

h(8) = subplot(3,3,8);
pl_L   = hist(bathy_data.h_L(:,2)/data.hm(2),bin_axis);
pl_rms = hist(bathy_data.h_krms(:,2)/data.hm(2),bin_axis);
ibeg = find(pl_L>0,1,'first'); iend = find(pl_L>0,1,'last');
plot( bin_axis(ibeg-1:iend+1)+df/2 , pl_L(ibeg-1:iend+1)/sum(pl_L(ibeg-1:iend+1)) , 'r' , 'linewidth' , 1 ); hold on, grid on, box on;
ibeg = find(pl_rms>0,1,'first'); iend = find(pl_rms>0,1,'last'); 
plot( bin_axis(ibeg-1:iend+1)+df/2 , pl_rms(ibeg-1:iend+1)/sum(pl_rms(ibeg-1:iend+1)) , 'k' , 'linewidth' , 1 )
% set(pl_rms,'edgecolor','k','linewidth',1)
% plot(data.hm,0.15,'s','markersize',4,'markerfacecolor',[0.8,0.5,0.1],'linewidth',1)
plot(median(bathy_data.h_krms(:,2)/data.hm(2)),0.015,'ko','markersize',4,'markerfacecolor','k')
pl_rms = errorbar(median(bathy_data.h_krms(:,2)/data.hm(2)),0.015,2*std(bathy_data.h_krms(:,2)/data.hm(2)),'horizontal','k','linewidth',0.75);
text(median(bathy_data.h_krms(:,2)/data.hm(2))-0.05,0.0098,'median','color','k','FontSize',8)
text(median(bathy_data.h_krms(:,2)/data.hm(2))-0.035,0.004,'\pm 2\sigma_h','color','k','FontSize',8)
plot(median(bathy_data.h_L(:,2)/data.hm(2)),0.015,'ro','markersize',4,'markerfacecolor','r')
pl_rms = errorbar(median(bathy_data.h_L(:,2)/data.hm(2)),0.015,2*std(bathy_data.h_L(:,2)/data.hm(2)),'horizontal','r','linewidth',0.75);
text(median(bathy_data.h_L(:,2)/data.hm(2))-0.06,0.0098,'median','color','r','FontSize',8)
text(median(bathy_data.h_L(:,2)/data.hm(2))-0.045,0.004,'\pm 2\sigma_h','color','r','FontSize',8)
set(gca,'Xtick',[0:0.1:4],'Xlim',[0.65 1.65],'Ylim',[0 0.06],'Ytick',[0:0.01:0.4]); set(gca, 'layer', 'top'), text(0.695,0.0525,'h)','Fontsize',10); hold off
% Statistics
text(0.82,0.0535,['Surveyed {\it h} = ',num2str(fix(data.hm(2)*100)/100),' m'],'color',[0.6,0.6,0.6],'FontSize',9)
text(0.687,0.042,['Median = ',num2str(fix(median(bathy_data.h_krms(:,2))*100)/100),' m'],'FontSize',8)
text(0.69,0.041-0.0065,['\sigma_h = ',num2str(fix(std(bathy_data.h_krms(:,2))*100)/100),' m'],'FontSize',8)
text(1.377,0.044,['Median = ',num2str(fix(median(bathy_data.h_L(:,2))*100)/100),' m'],'Color','r','FontSize',8)
text(1.38,0.043-0.0065,['\sigma_h = ',num2str(fix(std(bathy_data.h_L(:,2))*100)/100),' m'],'Color','r','FontSize',8)

h(9) = subplot(3,3,9);
pl_L   = hist(bathy_data.h_L(:,3)/data.hm(3),bin_axis);
pl_rms = hist(bathy_data.h_krms(:,3)/data.hm(3),bin_axis);
ibeg = find(pl_L>0,1,'first'); iend = find(pl_L>0,1,'last');
plot( bin_axis(ibeg-1:iend+1)+df/2 , pl_L(ibeg-1:iend+1)/sum(pl_L(ibeg-1:iend+1)) , 'r' , 'linewidth' , 1 ); hold on, grid on, box on;
ibeg = find(pl_rms>0,1,'first'); iend = find(pl_rms>0,1,'last'); 
plot( bin_axis(ibeg-1:iend+1)+df/2 , pl_rms(ibeg-1:iend+1)/sum(pl_rms(ibeg-1:iend+1)) , 'k' , 'linewidth' , 1 )
% set(pl_rms,'edgecolor','k','linewidth',1)
% plot(data.hm,0.15,'s','markersize',4,'markerfacecolor',[0.8,0.5,0.1],'linewidth',1)
plot(median(bathy_data.h_krms(:,3)/data.hm(3)),0.015,'ko','markersize',4,'markerfacecolor','k')
pl_rms = errorbar(median(bathy_data.h_krms(:,3)/data.hm(3)),0.015,2*std(bathy_data.h_krms(:,3)/data.hm(3)),'horizontal','k','linewidth',0.75);
text(median(bathy_data.h_krms(:,3)/data.hm(3))-0.05,0.0098,'median','color','k','FontSize',8)
text(median(bathy_data.h_krms(:,3)/data.hm(3))-0.035,0.004,'\pm 2\sigma_h','color','k','FontSize',8)
plot(median(bathy_data.h_L(:,3)/data.hm(3)),0.015,'ro','markersize',4,'markerfacecolor','r')
pl_rms = errorbar(median(bathy_data.h_L(:,3)/data.hm(3)),0.015,2*std(bathy_data.h_L(:,3)/data.hm(3)),'horizontal','r','linewidth',0.75);
text(median(bathy_data.h_L(:,3)/data.hm(3))-0.06,0.0098,'median','color','r','FontSize',8)
text(median(bathy_data.h_L(:,3)/data.hm(3))-0.045,0.004,'\pm 2\sigma_h','color','r','FontSize',8)
set(gca,'Xtick',[0:0.1:4],'Xlim',[0.65 1.65],'Ylim',[0 0.06],'Ytick',[0:0.01:0.4]); text(0.695,0.0525,'i)','Fontsize',10); hold off
% Statistics
text(0.82,0.0535,['Surveyed {\it h} = ',num2str(fix(data.hm(3)*100)/100),' m'],'color',[0.6,0.6,0.6],'FontSize',9)
text(0.687,0.042,['Median = ',num2str(fix(median(bathy_data.h_krms(:,3))*100)/100),' m'],'FontSize',8)
text(0.69,0.041-0.0065,['\sigma_h = ',num2str(fix(std(bathy_data.h_krms(:,3))*100)/100),' m'],'FontSize',8)
text(1.317,0.044,['Median = ',num2str(fix(median(bathy_data.h_L(:,3))*100)/100),' m'],'Color','r','FontSize',8)
text(1.32,0.043-0.0065,['\sigma_h = ',num2str(fix(std(bathy_data.h_L(:,3))*100)/100),' m'],'Color','r','FontSize',8)

% Axes
for hh = 7:9
  subplot(h(hh))
  set(gca, 'layer', 'top'); 
  set(gca,'TickDir','out','FontSize',9);
  xlabel('$h/h_{\rm survey}$ [-]','Interpreter','latex','Fontsize',12)
  xtickangle(0)
end

% Position and ticks
set(h(1),'Position',[0.06 0.728 0.27 0.26])
set(h(2),'Position',[0.40 0.728 0.27 0.26])
set(h(3),'Position',[0.72 0.728 0.27 0.26])
set(h(4),'Position',[0.06 0.39 0.27 0.28])
set(h(5),'Position',[0.40 0.39 0.27 0.28])
set(h(6),'Position',[0.72 0.39 0.27 0.28])
set(h(7),'Position',[0.06 0.09 0.27 0.22])
set(h(8),'Position',[0.40 0.09 0.27 0.22])
set(h(9),'Position',[0.72 0.09 0.27 0.22])

% Saving
% print(gcf,'-depsc','-r300',['inversion_illustration_1209_1900'])
print(gcf,'-dpng','-r300',['inversion_illustration_1209_1900.png'])
