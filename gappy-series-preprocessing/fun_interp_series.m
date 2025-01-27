function [ x , data ] = fun_interp_series( x , data , tail_method )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear interpolation of a timeseries, with special treatment at the timeseries tails.
%
% Inputs:
%   x    - grid data (generally time)
%   data - series of data containing NaNs
%   tail_method : 1 - imposing first (last) non-NaN value at beginning (end) of series
%                 2 - removing NaN values at beginning and end of series (NB: thus changing length of series)
%
% Outputs: 
%   x    - grid data, adjusted in case tail_method == 2
%   data - interpolated series
%
% June 27, 2023
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % ---------------------------------------------
  % Checking on inputs

  lx = length(x); x = x(:); data = data(:);
  if (lx ~= length(data)) error('Inputs have different sizes.'); end
  if (exist('tail_method') ~= 1)  tail_method = 1; end
  if (tail_method > 2) error('Unknown tail method.'); end
  

  % ---------------------------------------------
  % Dealing with beginning and ending sections of the series  
  
  idb = find( ~isnan(data) , 1 , 'first' ); % Beginning
  ide = find( ~isnan(data) , 1 , 'last' );  % End
  if (tail_method == 1)
    if (and(~isempty(idb),idb>1)), data(1:idb-1) = data(idb); end
    if (and(~isempty(ide),ide<lx)), data(ide+1:end) = data(ide); end
  else
    x = x(idb:ide); 
    data = data(idb:ide);
  end
  

  % ---------------------------------------------
  % Interpolation
  
  inNaNs = find( ~isnan(data) );
  if ~isempty(inNaNs)
    data   = interp1( x(inNaNs) , data(inNaNs) , x , 'linear' );
  end
  
  return
end
