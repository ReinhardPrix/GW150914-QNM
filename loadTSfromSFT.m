## Copyright (C) 2016 Reinhard Prix
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

function ts = loadTSfromSFT ( SFTpath, fMin, fMax, fSamp = [] )
  %% ts = loadTSfromSFT ( SFTpath, fMin, fMax, fSamp )
  %% load SFT with path 'SFTpath' and extract timeseries narrow-banded to [fMin,fMax]
  %% special values 'fMin/fMax = -1' mean to use up to the lowest/highest frequency of the SFT
  %% If 'fSamp' is not given, it defaults to the Nyquist rate fSamp = 2*fMax
  global debugLevel = 1;

  assert ( (fMin == -1) || (fMax == -1) || (fMax > fMin) );
  assert ( isempty ( fSamp ) || (fSamp >= 2*(fMax-fMin)), "Sampling frequency '%f' too low for bandwith '%f/2': will get aliasing\n", fSamp, 2*(fMax-fMin) );

  %% ----- read FT-data from SFT ----------
  sft = readSFT ( SFTpath );

  Tsft    = sft.header.Tsft;
  df      = 1/Tsft;
  Nfreq   = length ( sft.SFTdata );
  fMinSFT = sft.header.f0;
  fMaxSFT = fMinSFT + (Nfreq-1) * df;

  fkSFT   = fMinSFT + df * [ 0 : (Nfreq-1) ];
  xkSFT   = sft.SFTdata(:,1) + I * sft.SFTdata(:,2);

  ft0.fk    = fkSFT;
  ft0.xk    = xkSFT;
  ft0.IFO   = sft.header.IFO;
  ft0.epoch = sft.header.epoch.gpsSeconds;

  %% ----- handle special input values and defaults
  if ( fMin == -1 ) fMin = fMinSFT; endif
  if ( fMax == -1 ) fMax = fMaxSFT; endif
  if ( isempty ( fSamp ) ) fSamp = 2 * fMax; endif

  %% ----- turn FT into narrow-banded time-series at requested sampling rate
  ts = freqBand2TS ( ft0, fMin, fMax, fSamp );

  return;

endfunction %% loadTSfromSFT()
