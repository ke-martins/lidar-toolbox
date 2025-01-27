function [ s ] = fun_gridded_lidar_diagnostics( sf , x , data , bplot )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing basic statistics of gridded lidar data.
%
% Inputs:
%   sf    - sampling frequency [Hz]
%   x     - cross-shore grid [m]
%   data  - gridded data (dimensions: (t,x), typically the field 'z' output by fun_multibeam_lidar_gridding)
%   bplot - 1 for yes, 0 for no (default)
%
% Outputs: 
%   s     - a self-explanatory data structure containing some basic info on lidar data
%
% September 19, 2022
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Basic check on inputs
  if size(data,2) ~= length(x), error('Check inputs, the size of data array does not match the provided grid.'), end
  if nargin < 3, bplot = 0; end

  % Initialisation
  max_nb_points      = size(data,1); % maximal number of non-NaN points
  s.x                = x;
  s.x_info           = 'Original cross-shore grid [m]';
  s.per_returns      = nan(1,length(x));
  s.per_returns_info = 'Percentages of returns [%], i.e. % of non-NaN points';
  s.max_gap          = zeros(1,length(x))+1/sf;
  s.max_gap_info     = 'Largest gaps of data [s]';
  s.zeta_mean        = nanmean(data,1);
  s.zeta_mean_info   = 'Mean surface elevation [m]';
  % s.Hsig             = 4*nanstd(data-nanmean(data,1),[],1);
  s.Hsig             = 4*nanstd(data-nanmean(data,1),1);
  s.Hsig_info        = 'Significant wave height [m]';
  s.zeta_95p         = nan(1,length(x));
  s.zeta_95p_info    = 'Surface elevation 95% percentile [m] (top extreme values)';
  s.zeta_5p          = nan(1,length(x));
  s.zeta_5p_info     = 'Surface elevation 5% percentile [m] (top extreme values)';

  % First loop over spatial grid
  for xx = 1:length(x)
    % Computing percentages of returns
    s.per_returns(xx) = length(find(~isnan(data(:,xx))))/max_nb_points*100;

    % Computing the longest gap in data at this location
    % Note that it does not tell you much about spatial patterns
    iNaNs = isnan(data(:,xx));                       % Flagging NaNs
    iloc = reshape(find(diff([0;iNaNs;0])~=0),2,[]); % Reorganising by pair beginning and end of NaN zone
    [tmp,~] = max(diff(iloc));
    if ~isempty(tmp), s.max_gap(xx) = tmp/sf; end

    % Computing surface elevation distribution
    sorted_zeta    = sort(data(~isnan(data(:,xx)),xx));
    nb_points      = length(sorted_zeta);
    if nb_points > 2
      s.zeta_5p(xx)  = mean(sorted_zeta(1:fix(nb_points/20)));
      s.zeta_95p(xx) = mean(sorted_zeta(end-fix(nb_points/20):end));
    end
  end

  % Plot
  if bplot
    scrsz = get(0,'ScreenSize'); fig1 = figure;
    set(fig1,'Position',[50 50 scrsz(3)*0.5 scrsz(4)*0.55],'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 21 16],'color','w');
    set(0,'defaultAxesFontSize',9); set(0,'defaultaxeslinewidth',1)

    h(1) = subplot(2,2,1);
    plot( s.x, s.per_returns , 'k' , 'linewidth' , 1. ); grid on, box on;
    % Ticks
    ylabel('Returns [%]'), set(gca,'TickDir','out'); set(gca,'Layer','top')
    set(gca, 'xlim', [s.x(1) s.x(end)])
    set(gca, 'ytick', [0:20:100]), set(gca, 'ylim', [0 100])

    h(2) = subplot(2,2,2); hl = nan(1,5);
    plot( s.x, s.max_gap , 'k' , 'linewidth' , 1. ); grid on, box on;
    % Ticks
    ylabel('Maximal gap [s]','Interpreter','Latex','Fontsize',12), set(gca,'TickDir','out'); set(gca,'Layer','top')
    set(gca, 'xlim', [s.x(1) s.x(end)])

    h(3) = subplot(2,2,3);
    plot( s.x, s.Hsig , 'k' , 'linewidth' , 1. ); grid on, box on;
    % Ticks
    xlabel('$x$ [m]','Interpreter','Latex','Fontsize',12)
    ylabel('H_{sig} [m]'), set(gca,'TickDir','out'); set(gca,'Layer','top')
    set(gca, 'xlim', [s.x(1) s.x(end)])

    h(4) = subplot(2,2,4); hl = nan(1,3);
    hl(1) = plot( s.x , s.zeta_mean , 'color' , 'k' , 'linewidth' , 1. ); hold on, grid on, box on;
    hl(2) = plot( s.x , s.zeta_5p , '--' , 'color' , 'r' , 'linewidth', 1. );
    hl(3) = plot( s.x , s.zeta_95p , '--' , 'color' , [0.2,0.2,0.9] , 'linewidth', 1. );
    % Legend
    leg = legend(hl(1:3),'Mean','$\zeta_{5\%}$','$\zeta_{95\%}$','Location','NorthWest','Interpreter','Latex','Fontsize',10); legend boxoff
    leg.ItemTokenSize = [16,16,16]; set(leg,'FontSize',9)
    % Ticks
    xlabel('$x$ [m]','Interpreter','Latex','Fontsize',12)
    ylabel('$\zeta$ [m]','Interpreter','Latex','Fontsize',12), set(gca,'TickDir','out'); set(gca,'Layer','top')
    set(gca, 'xlim', [s.x(1) s.x(end)])

    % Positions
    set(h(1),'Position',[0.08 0.56 0.38 0.37])
    set(h(2),'Position',[0.54 0.56 0.38 0.37])
    set(h(3),'Position',[0.08 0.12 0.38 0.37])
    set(h(4),'Position',[0.54 0.12 0.38 0.37])
  end
end
