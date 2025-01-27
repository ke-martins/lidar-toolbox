function [ pNaNs ] = fun_count_pNaNs( data )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returning the percentage of NaNs in a timeseries.
%
% Inputs:
%   data - series of data (potentially) containing NaNs
%
% Outputs: 
%   pNaNs - percentage of NaNs in timeseries.
%
% June 27, 2023
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % ---------------------------------------------
  % Checking on inputs

  lx = length(data); data = data(:);
  isNaNs = find(isnan(data));
  pNaNs  = length(isNaNs)/lx*100;
  
  return
end
