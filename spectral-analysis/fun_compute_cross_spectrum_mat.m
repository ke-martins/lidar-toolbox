function [ data ] = fun_compute_cross_spectrum_mat( x , y , fs , wind )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of cross-spectrum [m^2/Hz] using a direct fft-based approach.
% The function is written for data provided in matrix (already organised).
% Overlapping (if any) thus has already been taken care of.
%
% Inputs:
%  x       - first signal  (organised in matrix length(timeseries),nb_block))
%  y       - second signal (organised in matrix length(timeseries),nb_block))
%  fs      - sampling frequency
%  wind    - Type of window for tappering ('rectangular', 'hann' or 'kaiser')
%
% Outputs: 
%  data    - a self-explanatory data structure containing spectra products
%             For more details, see through the code, where the data is stored.
%
% March 6, 2024
% Kévin Martins - kevin.martins@cnrs.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % --------------------- Various parameters -------------------------
  
  if (abs(size(x)/size(y)-1)>0.00000001), error('Wrong inputs.'); end
  if (isempty(wind) == 1) wind = 'rectangular'; end
  nfft   = size(x,1);
  nblock = size(x,2);


  % ---------------------- Initialization ------------------------
  % notes: nfft is forced to be even, this will be useful to get frequencies
  %        centered around 0.

  if (rem(nfft,2) ~= 0)
    x    = x(1:end-1,:);
    y    = y(1:end-1,:);
    nfft = nfft-1;
  end
  freqs = [-nfft/2:nfft/2]'/nfft*fs;

  % Output fields
  data.info      = 'Energy density spectra';
  data.f         = freqs;
  data.f_info    = 'Frequency [Hz] (one-sided)';
  data.df        = abs(freqs(2)-freqs(1));
  data.df_info   = 'Frequency resolution [Hz]';
  data.Exx       = zeros(nfft+1,1);
  data.Exx_info  = '(half) Energy density spectrum of first signal [m^2/Hz]';
  data.Eyy       = zeros(nfft+1,1);
  data.Eyy_info  = '(half) Energy density spectrum of second signal [m^2/Hz]';
  data.cpsd      = zeros(nfft+1,1);
  data.cpsd_info = 'Cross-spectrum between two signals [m^2/Hz]';
    
  % Local stuff
  Ax = zeros(nfft+1,nblock); % Fourier coefficients, for each block
  Ay = zeros(nfft+1,nblock); % Fourier coefficients, for each block
  nmid = (nfft)/2 + 1;       % Middle frequency (f = 0)


  % ---------------------- Compute FFT ----------------------------
  
  % Preparation of time windows
  switch wind
    case 'hann'
      ww = window(@hann,nfft); normFactor = mean(ww.^2);
    case 'kaiser'
      ww = window(@kaiser,nfft,3.5); normFactor = mean(ww.^2);
    case 'rectangular'
      ww = window(@rectwin,nfft); normFactor = mean(ww.^2);
  end

  % Computing FFT (loop over blocks)
  for kk = 1:nblock
    % Preparing block (window type)
    xseg = x(:,kk); xseg = (xseg(:) - mean(xseg)); % De-mean
    yseg = y(:,kk); yseg = (yseg(:) - mean(yseg)); % De-mean
    xseg = xseg.*ww / sqrt(normFactor);
    yseg = yseg.*ww / sqrt(normFactor);
    
    % FFT
    Ax_loc = fft( xseg , nfft )/nfft; Ay_loc = fft( yseg , nfft )/nfft;
    Ax(:,kk) = [ Ax_loc(nmid:nfft,:) ; Ax_loc(1:nmid,:) ]; % FFTshift
    Ay(:,kk) = [ Ay_loc(nmid:nfft,:) ; Ay_loc(1:nmid,:) ]; % FFTshift
    Ax(nmid,kk) = 0; Ay(nmid,kk) = 0;
  end

  % Accumulating products (loop over blocks)
  for kk = 1:nblock
    % Block kk FFT
    Ax_loc  = Ax(:,kk);
    Ay_loc  = Ay(:,kk);
    CAy_loc = conj(Ay(:,kk));
    
    % Compute PSDs and CPSD
    data.Exx = data.Exx + abs(Ax_loc.^2);
    data.Eyy = data.Eyy + abs(Ay_loc.^2);
    data.cpsd = data.cpsd + Ax_loc.*CAy_loc;
  end

  % Expected values
  data.Exx = data.Exx / nblock;
  data.Eyy = data.Eyy / nblock;
  data.cpsd = data.cpsd / nblock;
  
  
  % ---------------------- Finalisation ----------------------------

  % Keeping one-sided spectrum
  data.f = data.f(nmid:end);
  data.Exx  = data.Exx(nmid:end)/data.df;
  data.Eyy  = data.Eyy(nmid:end)/data.df;
  data.cpsd = data.cpsd(nmid:end)/data.df;
  
  % Number of blocks used to compute PSD
  data.nblocks      = nblock;
  data.nblocks_info = 'Number of blocks used to compute the PSD';

  return
end
