% Workflow example for the spectral analysis of gridded multibeam lidar data. 
% We start from gridded dataset, prepared in the "multibeam_lidar_gridding.m" file.
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

%% 1 - Spectral and bispectral analysis
% Re loading data files in case you haven't executed the code above
outfilename = 'data/BELS_20220912_1900_flight_2_x=225.dx=0.2m_2Hz_twin=0.16s.mat';
grid_data   = load(outfilename);

% Dealing with survey
survey = load('data/BELS_FRF_crawler+crab_survey_20220912_NAVD88.mat');

% Dealing with survey and water depth
interFunction = scatteredInterpolant([survey.x,survey.x]',[survey.y(1,:),survey.y(2,:)]',[survey.z(1,:),survey.z(2,:)]','linear','none');
grid_data.zb   = interFunction(grid_data.x,grid_data.y); clear interFunction
grid_data.zb_info = 'Interpolated seabed elevation [m] above NAVD88 datum'; 

% x = 220 m
[~,iloc] = nanmin(abs(grid_data.x-230));

% Spectral analysis
% Parameters for FFT and quality checks
sf = 2;
nfft = 128*sf;  % Number of points used per FFT - i.e. block length
overlap = 75;   % Amount of overlap in % between blocks of data
wind = 'hann';  % Tapering: use hann or rectangular
thperNaN = 10;  % 10 to 15% is completely acceptable from first tests

% Re-organising data by block
[~,zeta_mat] = fun_prep_gappy_series_by_block( grid_data.time, grid_data.z(:,iloc) , nfft , overlap , thperNaN );

% Computing PSD
psd_zeta = fun_compute_spectrum_mat( zeta_mat , sf , overlap , wind );

% Computing bispectrum
bis_zeta  = fun_compute_bispectrum_mat( zeta_mat , sf , overlap , 'rectangular' );

% Computing root-mean square wavenumber
h0 = nanmean(grid_data.z(:,iloc)) - grid_data.zb(iloc);
bis_zeta.k_rms      = fun_compute_krms( h0 , bis_zeta.f , bis_zeta.P , bis_zeta.B );
bis_zeta.k_rms_info = 'Boussinesq estimate of root-mean square wavenumber (Herbers et al., 2002)';

%% 2 - Plot of PSD and wave phase velocity spectrum
% Figure
scrsz = get(0,'ScreenSize'); fig1 = figure(1); 
set(fig1,'Position',[500 350 scrsz(3)*0.35 scrsz(4)*0.40],'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 20 10],'color','w');
set(0,'defaultAxesFontSize',8)

% K_L -- Approximation by Guo (2002) of the linear wave dispersion
nmid = (length(bis_zeta.f)-1)/2 + 1; % Middle frequency (f = 0)
kL   = (2*pi*bis_zeta.f(nmid:end)).^2/9.81 .* (1-exp(-((2*pi*bis_zeta.f(nmid:end))*sqrt(h0/9.81)).^(5/2))).^(-2/5);
 
% Wave phase speed spectra
h(1) = subplot(2,2,1); hl_1 = nan(1,3);
hl_1(1) = plot( bis_zeta.f , 0*bis_zeta.f + sqrt(9.81*h0) ,'k--','LineWidth',0.5); hold on, grid on, box on
hl_1(2) = plot( bis_zeta.f(nmid:end) , 2*pi*bis_zeta.f(nmid:end) ./ kL , 'r', 'LineWidth', 1);
hl_1(3) = plot( bis_zeta.f , 2*pi*bis_zeta.f ./ bis_zeta.k_rms , 'ko', 'markersize', 2., 'LineWidth', 0.5 ); hold off
set(gca, 'xlim', [0 0.5]), set(gca, 'xtick', [0 0.05 0.1:0.05:1],'Fontsize',9)
set(gca, 'ylim', [2 7]), set(gca, 'ytick', [0:1:30],'Fontsize',9)
ylabel( '$c(f) = 2\pi f/\kappa$ \,[m/s]', 'Interpreter', 'Latex', 'Fontsize', 11); 
set(gca,'TickDir','out');
text(0.012,6.6,'(a)','Fontsize',9,'FontWeight','bold')
xtickangle(0)
% Legend
leg = legend( hl_1 , '$\sqrt{gh_0}$', '$2\pi f/\kappa_L$','$2\pi f/\kappa_{rms}$','Location','South','Interpreter','Latex'); leg.ItemTokenSize = [16,16];
% leg = legend( hl_1 ,'$2\pi f/\kappa_L$','$2\pi f/\kappa_B$','$2\pi f/\kappa_{rms}$','Location','NorthEast'); leg.ItemTokenSize = [16,16,16];
leg_pos = get(leg,'Position'); leg_pos(1) = leg_pos(1)-0.02; leg_pos(2) = leg_pos(2)+0.04; set(leg,'Position',leg_pos);
set(leg,'Fontsize',10), legend boxoff

% PSD - A3
h(2) = subplot(2,2,3);
semilogy( psd_zeta.f, psd_zeta.E , 'k' , 'LineWidth' , 1 ); hold on, grid on, box on
set(gca, 'xlim', [0 0.5]), set(gca, 'xtick', [0 0.05 0.1:0.05:1],'Fontsize',9)
set(gca, 'ylim', [10^-3 10^1]), set(gca, 'ytick', [10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2],'Fontsize',9)
xlabel( '$f$ [Hz]', 'Interpreter', 'Latex', 'Fontsize', 11)
ylabel( '$E(f)$ [m$^2$/Hz]','Interpreter','Latex','Fontsize',11)
set(gca,'TickDir','out');
text(0.0105,4.25*10^0,'(b)','Fontsize',9,'FontWeight','bold')

% Add confidence levels 
coef_ulimit = psd_zeta.CI(2); coef_llimit = psd_zeta.CI(1);
f_CI = 0.075; psd_CI = 0.007;
std_upper_limit = psd_CI*coef_ulimit;
std_lower_limit = psd_CI*coef_llimit;
errorbar( f_CI , psd_CI , abs(psd_CI-std_lower_limit) , abs(psd_CI-std_upper_limit) , 'ko' , 'markerfacecolor' , 'k' , 'capsize' , 4 , 'markersize' , 2 , 'linewidth', 0.5)
text( 0.05 , 0.025 , '95% C.I.','Fontsize',7), pause(0.5)
xtickangle(0)

% Bicoherence (showing only a fraction, due to bispectrum symmetrical properties)
h(3) = subplot(2,2,2);
fid = find(and(bis_zeta.f>=0,bis_zeta.f<=0.5));
f_plot = bis_zeta.f(fid);
bic = bis_zeta.Bic(fid,fid);
ifc = find(f_plot>0.25,1,'first');
for ff = 2:ifc
  bic(ff+1:end,ff) = NaN;
end
for ff = ifc:2*ifc
  if ff <= numel(fid)
    bic(2*ifc-ff+1:end,ff) = NaN;
  end
end
[~,hcon,~] = contourf( f_plot , f_plot , bic.^2 ); hcon.LineWidth = 0.25; hold on, grid on
set(gca, 'xlim', [0 0.5]), set(gca, 'xtick', [0 0.1:0.1:1],'Fontsize',9)
set(gca, 'ylim', [0 0.275001]), set(gca, 'ytick', [0 0.05:0.05:1],'Fontsize',9)
xlabel( '$f$ [Hz]', 'Interpreter', 'Latex', 'Fontsize', 11)
ylabel( '$f$ [Hz]', 'Interpreter', 'Latex', 'Fontsize',11)
set(gca,'TickDir','out'); %set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');
text(0.015 , 0.262 ,'(c)','Fontsize',9,'FontWeight','bold')
x = [0.06 0.08];    % adjust length and location of arrow 
y = [0.13 0.09];      % adjust hieght and width of arrow
annotation('textarrow',[0.675 0.672]-0.03,...
  [0.415 0.3]+0.02,'String',' ($f_p,f_p$) ','FontSize',10,'Linewidth',0.75,'interpreter','latex','HeadWidth',4,'HeadLength',4)
xtickangle(0)

% Color bar
cmap = hot(128); colormap(flipud(cmap)); gca_pos = get(gca,'Position');
caxis([0 0.6]); pause(0.5)
hc = colorbar('XTick',0:0.1:0.6,'Location','northoutside');  
clabel = ylabel(hc,'$b^2(f,f)$ [-]','interpreter','latex','fontsize',10);
set(hc,'LineWidth',0.75); cl_pos = get(hc,'Position');
cl_pos(1) = cl_pos(1)+0.09;
cl_pos(2) = cl_pos(2)+0.03;
cl_pos(3) = 0.80*cl_pos(3);
cl_pos(4) = 0.5*cl_pos(4);
set(hc,'Position',cl_pos), pause(0.5)
set(gca,'Position',gca_pos)

% Positions
set(h(1),'Position',[0.08 0.60 0.41 0.375])
set(h(2),'Position',[0.08 0.14 0.41 0.375])
set(h(3),'Position',[0.595 0.14 0.395 0.69])

% Saving
print(fig1,'-dpng','-r300','PSD_krms_and_Bic_BELS_flight_#2.png')
