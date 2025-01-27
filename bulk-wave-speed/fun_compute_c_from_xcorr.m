function [ c , lag ] = fun_compute_c_from_xcorr( zeta_2 , zeta_1 , sf , dx , maxlag )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function computing 'bulk' wave propagation speed from cross-correlation of adjacent surface elevation timeseries.
% We use spline interpolation to reach sub-timestep resolution.
% For now, this assumes no NaNs in the provided timeseries.
%
% Inputs: 
%   zeta_2 - surface elevation timeseries at position 2 [m]
%   zeta_1 - surface elevation timeseries at position 1 [m]
%   sf     - sampling frequency [Hz]
%   dx     - distance between two locations where zeta is provided [m]
%   maxlag - max lag authorised for cross-correlation analysis [s]
%
% Outputs: 
%   c      - mean wave speed [m/s]
%   lag    - corresponding time lag [s]
%
% April, 2021
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --- Checking entries ---
  if (nargin == 4)
    warning('Default parameters are being used for the maxlag (10% the length of provided timeseries).')
    maxlag = fix(0.1*length(zeta_1)/sf);
  elseif or(nargin < 4,nargin > 5)
    error('Error: check entries.')
  end  
  iNaNs = [find(isnan(zeta_1)),find(isnan(zeta_2))];
  if (~isempty(iNaNs))
    c = NaN; lag = NaN;
    return
    warning('NaNs were detected in the timeseries, please deal with it first.')
  end

  % Cross-correlation
  [r,lags] = xcorr( zeta_2 , zeta_1 , fix(maxlag*sf) , 'coeff' );
  
  % Rough estimate
  [~,imaxcorr] = max(r);
  sub_lags = lags(imaxcorr-6:imaxcorr+6);
  sub_r    = r(imaxcorr-6:imaxcorr+6);
  
  % Looking for finer estimate using spline interpolation
  scale = 400;
  fine_lags = sub_lags(1):1/scale:sub_lags(end);
  fine_r    = interp1( sub_lags , sub_r , fine_lags , 'spline' );
  [~,imaxcorr] = max(fine_r);
  lag = (fine_lags(imaxcorr)/sf);
  c = dx/lag;
%   figure(2), plot( fine_lags , fine_r , 'k' ), hold on
end
