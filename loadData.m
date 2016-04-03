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

function [multiTS, multiPSD] = loadData ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"SFTpaths", "cell" },				%% cell-array [over detectors] containing SFT paths (1 SFT per IFO)
                        {"TimeRange", "real,strictpos,vector" },	%% time-window [tStart, tEnd] of analysis segment to return TS for
                        {"fSamp", "real,strictpos,scalar", [] },	%% sampling frequency for TS
                        {"assumeSqrtSX", "real,strictpos,vector", [] },	%% use this value (per IFO) for (white-noise) PSD instead of estimate
                        {"injectSqrtSX", "real,positive,vector", [] }	%% replace detector-data by Gaussian white noise of given PSD (can be 0)
                      );
  assert ( isvector ( uvar.TimeRange ) && (length(uvar.TimeRange) == 2) );

  numIFOs = length ( uvar.SFTpaths );
  if ( !isempty ( uvar.assumeSqrtSX ) )
    assert ( length ( uvar.assumeSqrtSX ) == numIFOs );
  endif
  if ( !isempty ( uvar.injectSqrtSX ) )
    assert ( length ( uvar.injectSqrtSX ) == numIFOs );
  endif

  multiPSD = cell ( 1, numIFOs );
  multiTS  = cell ( 1, numIFOs );
  for X = 1 : numIFOs

    %% ----- load time-series (TS) data from SFT
    ts0 = loadTSfromSFT ( uvar.SFTpaths{X}, -1, -1, uvar.fSamp );	%% load full-frequency-band TS from SFT
    dt = 1 / ts0.fSamp;

    %% ----- replace TS data by Gaussian noise if 'injectSqrtSn' given:
    if ( !isempty ( uvar.injectSqrtSX ) )
      sigma = uvar.injectSqrtSX(X) / sqrt(2 * dt);		%% convert PSD into Gaussian std-dev
      if ( sigma > 0 )
        ts0.xi = normrnd ( 0, sigma, size(ts0.xi) );	%% replace data by Gaussian noise
      else
        ts0.xi = zeros ( size(ts0.xi) );	%% normrnd() returns NANs for sigma=0 ...
      endif
    endif %% if injectSqrtSX

    %% ----- find time-samples corresponding to analysis segment 'TimeRange'
    samplesSeg = find ( (ts0.ti >= (uvar.TimeRange(1) - ts0.epoch)) & (ts0.ti <= (uvar.TimeRange(2) - ts0.epoch)) );

    %% ----- estimate PSD (*excluding* analysis segment 'TimeRange')
    multiPSD{X} = estimatePSD ( ts0, [samplesSeg(1), samplesSeg(end)] );

    if ( !isempty ( uvar.assumeSqrtSX ) )
      multiPSD{X}.Sn(:) = uvar.assumeSqrtSX(X)^2 ;	%% assume this PSD value instead of the measured estimate
    endif

    %% ----- extract TS of analysis segment 'TimeRange'
    tiStart = ts0.ti ( samplesSeg(1) );	%% exact start-time sample

    multiTS{X}        = ts0;
    multiTS{X}.ti     = ts0.ti ( samplesSeg ) - tiStart;	%% shift 'ti' to start from 0
    multiTS{X}.epoch += tiStart;
    multiTS{X}.xi     = ts0.xi ( samplesSeg );

  endfor %% X = 1 : numIFOs

  return;

endfunction  %% loadData()
