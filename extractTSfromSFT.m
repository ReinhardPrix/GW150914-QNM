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

function [ts, psd] = extractTSfromSFT ( varargin )
  global debugLevel = 1;
  global psd_version = 1;
  global cleanLines = false;

  uvar = parseOptions ( varargin,
                        {"SFTpath", "char,vector" },
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"Twindow", "real,strictpos,scalar", 10 },	%% time-window +- to extract around the event
                        {"plotSpectrum", "bool", false },
                        {"fSamp", "real,positive,scalar", 2000*2 },	%% sampling rate of output timeseries
                        {"injectionSources", "struct", []},
                        {"assumeSqrtSn", "real,strictpos,scalar", [] }
                      );
  assert ( uvar.fMax > uvar.fMin );

  pieces = strsplit ( uvar.SFTpath, "/" );
  sftfname = pieces { end };
  sftbname = strrep ( sftfname, ".sft", "");
  bname = sprintf ( "%s-freq%.0fHz-%.0fHz-fSamp%.0fHz-GPS%.0fs+-%.0fs-psd_v%d-lineCleaning%s",
                    sftbname, uvar.fMin, uvar.fMax, uvar.fSamp, uvar.tCenter, uvar.Twindow, psd_version, ifelse ( cleanLines, "On", "Off" ) );

  resDir = "extractTS-Results";
  [status, msg] = mkdir ( resDir );
  assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );

  psd_fname = sprintf ( "%s/%s.psd", resDir, bname );
  ts_fname  = sprintf ( "%s/%s.ts",  resDir, bname );
  %% extract IFO name from SFT name: only works for SFT-name compliant SFTs
  pieces = strsplit ( sftfname, {"-", "_"} );
  IFO = pieces{3};
  assert ( (length ( IFO ) == 2) && isalpha(IFO(1)) && isdigit(IFO(2)) );

  %% ---------- otherwise: Read SFT frequency-domain data ----------
  DebugPrintf (2, "%s: Extracting TS from SFT '%s'\n", funcName(), uvar.SFTpath );
  sft = readSFT ( uvar.SFTpath );
  epoch = sft.header.epoch.gpsSeconds;
  f0 = sft.header.f0;
  assert ( strcmp ( IFO, sft.header.IFO ) );
  Tsft = sft.header.Tsft;
  df = 1/Tsft;
  Nfreq = length ( sft.SFTdata );
  f1 = f0 + (Nfreq-1) * df;
  fk0 = f0 + df * [ 0 : (Nfreq-1) ];
  xk0 = sft.SFTdata(:,1) + I * sft.SFTdata(:,2);
  assert ( length(fk0) == length(xk0) );
  ft0.fk = fk0;
  ft0.xk = xk0;
  ft0.IFO = IFO;
  ft0.epoch = epoch;

  %% ---------- compute PSD on short timeseries, nuke lines, extract 'physical' frequency band, and whiten + overwhitened TS ----------
  switch ( psd_version )
    case 1
      %% ---------- extract frequency band of interest [fMin,fMax] as a timeseries ----------
      tsBand0 = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, uvar.fSamp );
      %%figure(); clf; plot ( tsBand0.ti + (epoch - uvar.tCenter), tsBand0.xi, "-" );
      indsTrunc = find ( (tsBand0.ti >= (uvar.tCenter - tsBand0.epoch - uvar.Twindow)) & (tsBand0.ti <= (uvar.tCenter - tsBand0.epoch + uvar.Twindow )) );
      %% ---------- truncate timeseries to [ tCenter - dT, tCenter + dT ] ----------
      tsBand.ti    = tsBand0.ti ( indsTrunc );
      tsBand.xi    = tsBand0.xi ( indsTrunc );
      tsBand.IFO   = IFO;
      tsBand.epoch = epoch;
      [ts, psd] = whitenTS ( "tsIn", tsBand,
                             "fMin", uvar.fMin, "fMax", uvar.fMax,
                             "plotSpectrum", uvar.plotSpectrum );
      assert ( isempty ( uvar.injectionSources ), "Sorry, QNM injections only supported with 'psd_version=2'\n");
      assert ( isempty ( uvar.assumeSqrtSn ), "Sorry, 'signal-only' simulation only supported with 'psd_version=2'\n");

    case 2
      [ts, psd] = whitenTS_v2 ( "ftIn", ft0, ...
                                "fSamp", uvar.fSamp, ...
                                "tCenter", uvar.tCenter, "Twindow", uvar.Twindow, ...
                                "fMin", uvar.fMin, "fMax", uvar.fMax, ...
                                "plotSpectrum", uvar.plotSpectrum, ...
                                "injectionSources", uvar.injectionSources, ...
                                "assumeSqrtSn", uvar.assumeSqrtSn
                              );
    otherwise
      error ("psd_version = %d not supported\n", psd_version );
  endswitch

  if ( uvar.plotSpectrum )
    fname = sprintf ( "%s/%s-spectrum.pdf", resDir, bname);
    ezprint ( fname, "width", 512 );
  endif

  %% ---------- store results for potential future re-use ----------
  fid = fopen ( psd_fname, "wb" ); assert ( fid != -1, "Failed to open '%s' for writing\n", psd_fname );
  fprintf ( fid, "%%%%%14s %16s\n", "freq [Hz]", "SX [1/Hz]" );
  fprintf ( fid, "%16.9f  %16.9g\n", [psd.fk', psd.Sn']' );
  fclose(fid);

  fid = fopen ( ts_fname, "wb" ); assert ( fid != -1, "Failed to open '%s' for writing\n", ts_fname );
  fprintf ( fid, "%%%%%14s %16s %16s %16s %16s\n", "ti [offs s]", "xi", "xi/sqrtSX", "xi/SX", "epoch [GPS s]" );
  epochV = ts.epoch * ones ( size ( ts.ti ) );	%% stupid way, but ensure we keep epoch separate for 'offsets' ti to avoid roundoff problems
  fprintf ( fid, "%16.9f %16.9g %16.9g %16.9g %16.9f\n", [ts.ti', ts.xi', ts.xiW', ts.xiOW', epochV']' );
  fclose(fid);

  return;

endfunction %% extractTSfromSFT()
