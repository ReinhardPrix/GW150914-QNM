#!/usr/bin/octave -q

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


global debugLevel = 1;

%% --------------------------------------------------
%% run a ringdown search as function of start-time 't0' on GW150914
%% --------------------------------------------------

%% ========== driver parameters ==========
SFTs = {"./Data/H-1_H1_1800SFT_ER8-C01-1126257832-1800.sft"; "./Data/L-1_L1_1800SFT_ER8-C01-1126258841-1800.sft" };
fSamp = 4000;	%% full sampling frequency of fmax=2kHz SFT, and conveniently such that 7.0ms timeshift between IFOs
                %% can be represented by exactly by an integer bin-shift: 7e-3 s * 4e3 Hz =  28.0 bins
shiftL1 = 7.0e-3;	%% time-shift to apply to L1 data-stream: currently 'official' value (v8)

confidence = 0.90;
useTSBuffer = true;

doPlotSnapshots = false;
doPlotContours  = false;
doPlotSummary   = false;
doPlotSpectra   = false;
doPlotBSGHist   = false;
doPlotH         = false;

if ( !exist("searchType") )     searchType = "verify"; endif
if ( !exist("extraLabel") )     extraLabel = ""; endif
if ( !exist("psd_version") )    global psd_version = 2; endif
if ( !exist("cleanLines") )     global cleanLines = false; endif
if ( !exist("data_FreqRange") ) data_FreqRange  = [ 30, 1e3 ]; endif
if ( !exist("iFig0") )          global iFig0 = 0; endif
if ( !exist("injectionSources") ) injectionSources = []; endif

%% ----- 'GR predictions/values on GW150914 ----------
tMergerGW150914 = 1126259462.42285;	%% from https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TestingGR/O1/G184098/ringdown_presence
f0GR = struct ( "val", 251, "lerr", 9, "uerr", 5 );	%% taken from 'testing GR paper' v21: DCC https://dcc.ligo.org/LIGO-P1500213-v21
taumsGR = struct ( "val", 4, "lerr", 0.2, "uerr", 0.4 );
plotMarkers = struct ( "name", "IMR", "f0", f0GR.val, "tau", taumsGR.val * 1e-3 );	%% by default: show IMR parameters on PE plots
%% ---------- Prior range defaults ----------

prior_f0Range   = [ 200, 300 ];
prior_tauRange  = [ 0.5e-3, 20e-3 ];
if ( !exist ( "prior_H" ) )
  H_i = [ 2e-22, 3e-22, 4e-22, 5e-22, 6e-22, 7e-22, 8e-22, 9e-22, 10e-22 ];
  p_i = 1 ./ H_i;
  p_i /= sum ( p_i(:) );
  prior_H = [ H_i', p_i' ];
endif
step_f0         = 0.5;
step_tau        = 0.5e-3;

%% NR best-fit 'merger time'
%% from https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TestingGR/O1/G184098/ringdown_presence

switch ( searchType )
  case "verify"
    %% ---------- test-case to compare different code-versions on ----------
    tMerger = tMergerGW150914;
    tOffsV = [7e-3];

    useTSBuffer = false;

    doPlotContours  = true;
    doPlotSummary   = false;
    doPlotSnapshots = true;
    doPlotSpectra   = true;
    doPlotH         = true;

  case "onSource"
    %% ---------- "ON-SOURCE ----------
    tMerger = tMergerGW150914;
    tOffsV = [ 1, 3, 5, 7 ] * 1e-3;

    useTSBuffer = false;

    doPlotSnapshots = true;
    doPlotContours  = true;
    doPlotSummary   = true;
    doPlotSpectra   = false;
    doPlotBSGHist   = false;
    doPlotH         = true;

  case "offSource"
    %% ---------- "OFF-SOURCE" for background estimation ----------
    tMerger = tMergerGW150914 + 10;
    tOffsV = [ -3 : 0.05 : 3 ];
    plotMarkers = [];

    useTSBuffer = true;

    doPlotSnapshots = true;
    doPlotContours  = false;
    doPlotSummary   = true;
    doPlotSpectra   = false;
    doPlotBSGHist   = true;

  case "injection"
    %% ---------- add QNM signal to test parameter-estimation accuracy ----------
    tMerger = tMergerGW150914 + 10;;	%% go to some 'off source' time-stretch
    tOffsV = [10e-3];
    injectionSources = struct ( "name", "Inj", "t0", tMerger + 10e-3, "A", 5e-22, "phi0", 0.5, "f0", 251, "tau", 4e-3, "shiftL1", shiftL1 );
    plotMarkers = injectionSources;

    useTSBuffer     = false;

    doPlotContours  = true;
    doPlotSummary   = false;
    doPlotSnapshots = true;
    doPlotSpectra   = false;
    doPlotH         = false;

  otherwise
    error ("Unknown searchType = '%s' specified\n", searchType );

endswitch

%% ----- data preparation -----

%% load frequency-domain data from SFTs:
for X = 1:length(SFTs)
  [ts{X}, psd{X}] = extractTSfromSFT ( "SFTpath", SFTs{X}, "fMin", min(data_FreqRange), "fMax", max(data_FreqRange), "fSamp", fSamp, ...
                                       "tCenter", fix(tMerger), "Twindow", 4, ...
                                       "plotSpectrum", doPlotSpectra, "useBuffer", useTSBuffer,
                                       "injectionSources", injectionSources
                                     );
  %% for plotting OW-timeseries: store noise-values at 'GR frequency f0GR'
  [val, freqInd] = min ( abs ( psd{X}.fk - f0GR.val ) );
  ts{X}.SX_GR = psd{X}.Sn ( freqInd );
endfor




%% create unique time-tagged 'ResultsDir' for each run:
gm = gmtime ( time () );
resDir = sprintf ( "Results/Results-%02d%02d%02d-%02dh%02d-%s-data%.0fHz-%.0fHz-psd_v%d-lineCleaning%s%s",
                   gm.year - 100, gm.mon + 1, gm.mday, gm.hour, gm.min, searchType, data_FreqRange,
                   psd_version, ifelse ( cleanLines, "On", "Off" ),
                   extraLabel );
[status, msg, id] = mkdir ( resDir ); assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );
addpath ( pwd() );
cd ( resDir );

%%try
%% ----- run search
Nsteps = length(tOffsV);

ret = cell ( 1, Nsteps );
for i = 1:Nsteps
  DebugPrintf ( 1, "t0GPS = tMerger + tOffs = %.6f s + %.1f ms\n", tMerger, tOffsV(i) * 1e3 );

  ret{i} = searchRingdown ( "ts", ts, "psd", psd, ...
                            "t0GPS", tMerger + tOffsV(i), ...
                            "prior_f0Range", prior_f0Range, "step_f0", step_f0, ...
                            "prior_tauRange", prior_tauRange, "step_tau", step_tau, ...
                            "prior_H", prior_H, ...
                            "shiftL1", shiftL1
                          );
  ret{i}.tMerger = tMerger;
  ret{i}.tOffs   = tOffsV(i);

  %% ----- save posterior in matrix format ----------
  fname = sprintf ( "%s-BSG.dat", ret{i}.bname );
  tmp = ret{i}.BSG_f0_tau;
  save ( "-ascii", fname, "tmp" );

  fname = sprintf ( "%s-SNR.dat", ret{i}.bname );
  tmp = ret{i}.SNR;
  save ( "-ascii", fname, "tmp" );

  %% ----- compute derived quantities
  %% 2D and 1D marginalized posteriors on {f,tau} ----------
  tmp = ret{i}.BSG_f0_tau;
  ret{i}.posterior2D = tmp / sum ( tmp(:) );

  tmp = sum ( ret{i}.posterior2D, 1 );
  ret{i}.posterior_f0   = tmp / sum ( tmp(:) );

  tmp = sum ( ret{i}.posterior2D, 2 );
  ret{i}.posterior_tau   = tmp / sum ( tmp(:) );

  f0  = unique ( ret{i}.ff0 );
  tau = unique ( ret{i}.ttau );
  ret{i}.f0_est   = credibleInterval ( f0,   ret{i}.posterior_f0,   confidence );
  ret{i}.tau_est  = credibleInterval ( tau,  ret{i}.posterior_tau,  confidence );
  ret{i}.isoConf2 = credibleContourLevel ( ret{i}.posterior2D, confidence );

  %% MPE values
  ret{i}.posteriorMax = max ( ret{i}.posterior2D(:) );
  ret{i}.l_MPE2    = l_MPE2 = ( find ( ret{i}.posterior2D(:) == ret{i}.posteriorMax ) )(1);
  ret{i}.f0_MPE2   = ret{i}.ff0 ( l_MPE2 );
  ret{i}.tau_MPE2  = ret{i}.ttau ( l_MPE2 );

  ret{i}.A_MPE2    = ret{i}.A_est( l_MPE2 );
  ret{i}.phi0_MPE2 = ret{i}.phi0_est ( l_MPE2 );
  ret{i}.SNR_MPE2  = ret{i}.SNR ( l_MPE2 );
  %% ---------- plot results summary page
  %%fname = sprintf ( "%s.png", ret{i}.bname );
  %%ezprint ( fname, "width", 1024, "height", 786, "dpi", 72 );
  if ( doPlotSnapshots )
    plotSnapshot ( ret{i}, [], plotMarkers );
  endif

  DebugPrintSummary ( 1, ret{i} );

endfor %% i = 1:Nsteps

%% ----- save axis {f0, tau} in matrix format ----------
tmp = ret{i}.ff0;
save ( "-ascii", "f0s.dat", "tmp" );
tmp = ret{i}.ttau;
save ( "-ascii", "taus.dat", "tmp" );

%% ----- store complete results dump ----------
fname = sprintf ( "RingdownDriver-%s.hd5", ret{1}.bname );
save ("-hdf5", fname )

%% ----- plot quantities vs tOffs ----------
if ( doPlotSummary )
  plotSummary ( ret );
endif

if ( doPlotContours )
  plotContours ( ret, [], plotMarkers )
endif

if ( doPlotH )
  figure(); clf; hold on;
  plot ( prior_H(:,1), prior_H(:,2), "-xg;prior;", "linestyle", "--", "linewidth", 2 );
  for i = 1 : Nsteps
    plot ( prior_H(:,1), ret{i}.post_H, "-b" );
    leg = sprintf ( "+%.1fms", ret{i}.tOffs * 1e3 );
    text ( prior_H(end,1), ret{i}.post_H(end), leg );
  endfor
  xlabel ( "H" ); ylabel ("pdf");
  grid on; hold off;
  fname = sprintf ( "%s-post_H.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );
endif

if ( doPlotBSGHist )
  figure ( iFig0  + 6 ); clf;
  BSG = zeros ( 1, Nsteps );
  for i = 1 : Nsteps
    BSG(i) = ret{i}.BSG;
  endfor
  hist ( BSG, 20 );
  xlabel ( "<BSG>" );
  fname = sprintf ( "%s-hist.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );
endif

cd ("../..");

## catch
##   err = lasterror();
##   warning(err.identifier, err.message);
##   err.stack.file
##   err.stack.name
##   err.stack.line

##   cd ("../..");	%% make sure we end up in main dir in case of failure
## end_try_catch


