## Copyright (C) 2015 Reinhard Prix
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

%% estimate single-sided noise PSD 'psd(f)', and return whitened and over-whitened TS
%% including an automated line-nuking algorithm: 'lines' are identified as > lineSigma deviations in power
function [ tsOut, psd ] = whitenTS ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"tsIn", "struct" },
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"lineSigma", "real,positive,scalar", 5},	%% sigma deviations to indentify 'lines' in spectrum
                        {"lineWidth", "real,positive,scalar", 0.1},	%% +- width in Hz to zero around 'lines'
                        {"RngMedWindow", "real,positive,scalar", 300 },  %% window size to use for rngmed-based PSD estimation
                        {"plotSpectrum", "bool", false }
                      );

  %% handy shortcuts
  epoch = uvar.tsIn.epoch;
  ti = uvar.tsIn.ti;
  xi = uvar.tsIn.xi;
  IFO = uvar.tsIn.IFO;

  %% re-normalize TS: make sure start-time = 0, fully carried by 'epoch'
  tStart = min ( ti );
  ti -= tStart;
  epoch += tStart;
  dt = mean ( diff ( ti ) );
  fSamp = round ( 1/dt );
  T = max(ti) + dt;
  ft = FourierTransform ( ti, xi );
  ft.IFO = IFO;
  ft.epoch = epoch;

  indsUse = binRange ( uvar.fMin, uvar.fMax, ft.fk );

  %% ---------- estimate PSD using running-median over frequency ----------
  %% Wiener-Khintchine theorm
  %% Sn_double = 1/T * mean ( periodo );
  if (uvar.RngMedWindow == 300 )
    rngmedbias = 0.694816624988959; %% XLALRngMedBias (300);	%% avoid needing swig-wrapped LAL for default
  else
    lal;
    rngmedbias = XLALRngMedBias ( uvar.RngMedWindow );
  endif

  periodo = abs ( ft.xk(indsUse) ).^2;
  Sn_wide = 2/T * rngmed ( periodo, uvar.RngMedWindow ) / rngmedbias;	%% single-sided PSD

  psd.fk = ft.fk ( indsUse );
  psd.Sn = Sn_wide;
  psd.IFO   = IFO;
  psd.epoch = epoch;

  %% ---------- automatically identify and nuke 'lines' ----------
  fk_wide   = ft.fk ( indsUse );
  xk_wide   = ft.xk ( indsUse );

  yk = xk_wide  ./ sqrt( T * Sn_wide / 2 );   %% normalized SFT:
  Pk_re = abs ( real ( yk ) );
  Pk_im = abs ( imag ( yk ) );
  inds_lines = find ( (Pk_re > uvar.lineSigma) | (Pk_im > uvar.lineSigma) );
  iNukeWidth = round ( uvar.lineWidth * T );
  sides = -iNukeWidth : iNukeWidth;
  [xx, yy] = meshgrid ( inds_lines, sides );
  nuke_lines = xx + yy;
  %%nuke_bands = find ( (fk_wide < 100) | (fk_wide >= 300 ) );
  %%nuke = [ nuke_lines(:); nuke_bands(:) ];
  nuke = [ nuke_lines(:) ];
  indsNuke = unique ( nuke(:) );
  indsNuke = indsNuke ( (indsNuke >= 1) & (indsNuke <= length(indsUse)) );

  %% ----- replace all lines with zeros ----------
  DebugPrintf ( 2, "\n----- Identified lines in %s: -----\n", IFO );
  DebugPrintf ( 2, "Line-frequencies: %f Hz\n", fk_wide ( inds_lines ) );

  %% ----- whitened and overwhitened sFTs (with nuked lines) ----------
  if ( uvar.plotSpectrum )
    if ( strcmp ( IFO, "H1" ) )
      psd0 = load ( "./Data/lalinferencemcmc-0-H1L1-1126259462.39-0H1-PSD.dat" );
    elseif ( strcmp ( IFO, "L1" ) )
      psd0 = load ( "./Data/lalinferencemcmc-0-H1L1-1126259462.39-0L1-PSD.dat");
    endif
    iIFO = ifelse ( strcmp ( IFO, "H1" ), 0, 1 );
    figure (); clf;

    subplot ( 3, 1, 1); hold on;
    plot ( fk_wide, abs ( xk_wide ) / sqrt(T), "+-", "color", "blue" ); legend ( IFO );
    if ( exist ( "psd0" ) )
      plot ( psd0(:,1), sqrt(psd0(:,2)), "x;lalinference;", "color", "magenta" );
    endif
    plot ( fk_wide, sqrt(Sn_wide), "o;sqrt(SX);", "color", "green" );
    plot ( fk_wide ( indsNuke ), abs ( xk_wide ( indsNuke )) / sqrt(T), "o", "color", "red" );
  endif

  NindsNuke = length ( indsNuke );
  draws = zeros ( 2, NindsNuke );
  Sn0 = Sn_wide ( indsNuke );
  noise = sqrt(T/4 * Sn0 ) .* ( draws(1,:) + I * draws(2,:) );
  xk_wide ( indsNuke ) = noise;
  xkW_wide  = xk_wide ./ sqrt(Sn_wide);
  xkOW_wide = xk_wide ./ Sn_wide;

  if ( uvar.plotSpectrum )
    plot ( fk_wide ( indsNuke ), abs ( xk_wide ( indsNuke )) / sqrt(T), "x", "color", "red" );
    xlim ( [uvar.fMin, uvar.fMax] );
    ylim ( [0, 5e-23] );
    grid on;
    hold off;

    subplot ( 3, 1, 2 );
    plot ( fk_wide, abs ( xkW_wide ) / sqrt(T/2), "+-", "color", "blue" ); legend ( "xk/rms" );
    xlim ( [uvar.fMin, uvar.fMax] );
    ylim ( [ 0, 5 ] );

    subplot ( 3, 1, 3 );
    plot ( fk_wide, abs ( xkOW_wide ), "+-", "color", "blue" ); legend ( "xk/Sn" );
    xlim ( [uvar.fMin, uvar.fMax] );
  endif

  %% ----- turn SFTs back into timeseries ----------
  ftOut.fk = fk_wide;
  ftOut.IFO = IFO;
  ftOut.epoch = epoch;
  ftOut.xk   = xk_wide;
  ftOut.xkW  = xkW_wide;
  ftOut.xkOW = xkOW_wide;

  ft0 = ftOut;
  ts   = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, fSamp );
  ft0.xk = ftOut.xkW;
  tsW  = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, fSamp );
  ft0.xk = ftOut.xkOW;
  tsOW = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, fSamp );

  assert ( (ts.ti == tsW.ti) && (ts.ti == tsOW.ti) );

  tsOut.ti   = ts.ti;
  tsOut.xi   = ts.xi;
  tsOut.xiW  = tsW.xi;
  tsOut.xiOW = tsOW.xi;
  tsOut.IFO  = IFO;
  tsOut.epoch= epoch;

  assert ( size(tsOut.ti) == size( ti ) );
  return;
endfunction
