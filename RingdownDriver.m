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
numIFOs = length ( SFTs );
fSamp = 4000;	%% full sampling frequency of fmax=2kHz SFT, and conveniently such that 7.0ms timeshift between IFOs
                %% can be represented by exactly by an integer bin-shift: 7e-3 s * 4e3 Hz =  28.0 bins
FreqRange  = [ 30, 1900 ]; %% SFT data spans [10,2e3]Hz, but avoid PSD 'ramp up/down' at boundaries

%% collect time-delay and antenna-pattern factors into one object
skyCorr = struct ( "IFO", {"H1", "L1"}, "timeShift", { 0, 7e-3 }, "ampFact", { 1, -1 } );	%% values specific to GW150914

confidence = 0.90;

doPlotSnapshots = false;
doPlotContours  = false;
doPlotSummary   = false;
doPlotSpectra   = false;
doPlotBSGHist   = false;
doPlotH         = false;
doPlotPErecovery= false;

if ( !exist("searchType") )     searchType = "verify"; endif
if ( !exist("extraLabel") )     extraLabel = ""; endif
if ( !exist("Tseg") )		Tseg = 8; endif				%% length of data segment (in s) to analyse, centered on 'tMerger'
if ( !exist("injectionSources") ) injectionSources = []; endif
if ( !exist("assumeSqrtSX") ) 	assumeSqrtSX = []; endif
if ( !exist("injectSqrtSX") ) 	injectSqrtSX = []; endif
if ( !exist("noMismatchInj") )  noMismatchInj = false; endif		%% inject at exact {f0,tau} grid-points for zero mismatch

%% ----- 'GR predictions/values on GW150914 ----------
tMergerGW150914 = 1126259462.42285;	%% from https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TestingGR/O1/G184098/ringdown_presence
f0GR = struct ( "val", 251, "lerr", 9, "uerr", 5 );	%% taken from 'testing GR paper' v21: DCC https://dcc.ligo.org/LIGO-P1500213-v21
taumsGR = struct ( "val", 4, "lerr", 0.2, "uerr", 0.4 );
plotMarkers = struct ( "name", "IMR", "f0", f0GR.val, "tau", taumsGR.val * 1e-3 );	%% by default: show IMR parameters on PE plots
%% ---------- Prior range defaults ----------

prior_f0Range   = [ 200, 300 ];
prior_tauRange  = [ 0.5e-3, 20e-3 ];
if ( !exist ( "prior_H" ) )
  H_i = [ 2 : 10 ] * 1e-22;
  p_i = 1 ./ H_i;
  p_i /= sum ( p_i(:) );
  prior_H = struct ( "x", H_i, "px", p_i );
endif
step_f0         = 0.5;
step_tau        = 0.5e-3;

%% NR best-fit 'merger time'
%% from https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TestingGR/O1/G184098/ringdown_presence

switch ( searchType )
  case "verify"
    %% ---------- test-case to compare different code-versions on ----------
    tMerger = tMergerGW150914;
    t0V = tMerger + [7] * 1e-3;

    doPlotContours  = false;
    doPlotSummary   = false;
    doPlotSnapshots = true;
    doPlotSpectra   = false;
    doPlotH         = false;

  case "onSource"
    %% ---------- "ON-SOURCE ----------
    tMerger = tMergerGW150914;
    t0V = tMerger + [ 1, 3, 5, 7 ] * 1e-3;

    doPlotSnapshots = true;
    doPlotContours  = true;
    doPlotSummary   = true;
    doPlotSpectra   = false;
    doPlotBSGHist   = false;
    doPlotH         = true;

  case "offSource"
    %% ---------- "OFF-SOURCE" for background estimation ----------
    tMerger = tMergerGW150914 + 10;
    numTrials = 100;
    t0V = tMerger + unifrnd ( -3, 3, 1, numTrials );
    plotMarkers = [];

    doPlotSnapshots = false;
    doPlotContours  = false;
    doPlotSummary   = true;
    doPlotSpectra   = false;
    doPlotBSGHist   = true;
    doPlotH         = false;

  case "injections"
    %% NOTE: injections defaults to injections into real-data, but override with 'injectSqrtSX' and 'assumeSqrtSX' settings
    tMerger = tMergerGW150914 + 10;	%% go to some 'off source' time-stretch
    numTrials = 10;
    t0V = tMerger + unifrnd ( -3, 3, 1, numTrials );
    injectionSources = struct();
    for m = 1 : numTrials
      injectionSources(m) = struct ( "name", 	sprintf("QNM-%d", m), ...
                                     "t0", 	t0V(m), ...
                                     "A", 	unifrnd ( 3e-22, 8e-22 ), ...
                                     "phi0", 	unifrnd ( 0, 2*pi ), ...
                                     "f0",	unifrnd ( min(prior_f0Range), max(prior_f0Range) ), ...
                                     "tau", 	unifrnd ( min(prior_tauRange), max(prior_tauRange) ), ...
                                     "skyCorr", skyCorr ...
                                   );
      if ( noMismatchInj )	%% place signal-{f0,tau} on exact search grid point
        injectionSources.f0  = min(prior_f0Range) + step_f0 * round ( (injectionSource.f0 - min(prior_f0Range)) / step_f0 );
        injectionSources.tau = min(prior_tauRange) + step_tau * round ( (injectionSource.tau - min(prior_tauRange)) / step_tau);
      endif
    endfor

    doPlotContours  = false;
    doPlotSummary   = false;
    doPlotSnapshots = false;
    doPlotSpectra   = false;
    doPlotH         = false;
    doPlotPErecovery= true;
    PErecovery = struct();

  otherwise
    error ("Unknown searchType = '%s' specified\n", searchType );

endswitch %% searchType


%% ========== create unique time-tagged 'ResultsDir' for each run:
gm = gmtime ( time () );
resDir = sprintf ( "Results/Results-%02d%02d%02d-%02dh%02d-%s-data%.0fHz-%.0fHz%s",
                   gm.year - 100, gm.mon + 1, gm.mday, gm.hour, gm.min, searchType, FreqRange, extraLabel );
if ( noMismatchInj )
  resDir = strcat ( resDir, "-noMismatchInj" );
endif
if ( !isempty ( injectSqrtSX ) )
  resDir = strcat ( resDir, "-injectSqrtSX", sprintf("_%.2g", injectSqrtSX ) );
endif
if ( !isempty ( assumeSqrtSX ) )
  resDir = strcat ( resDir, "-assumeSqrtSX", sprintf("_%.2g", assumeSqrtSX ) );
endif
[status, msg, id] = mkdir ( resDir ); assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );

bname = sprintf ( "Ringdown-f%.0fHz-%.0fHz-tau%.1fms-%.1fms",
                  min(prior_f0Range), max(prior_f0Range),
                  1e3 * min(prior_tauRange), 1e3 * max(prior_tauRange)
                );

%% ========== start actual analysis ====================

%% ----- data reading + preparation -----
TimeRange = [ tMerger - 0.5 * Tseg, tMerger + 0.5 * Tseg ];
DebugPrintf ( 1, "Loading SFT data + PSD-estimation ... ");
[multiTS0, multiPSD0] = loadData ( "SFTpaths", SFTs, "TimeRange", TimeRange, "fSamp", fSamp, "assumeSqrtSX", assumeSqrtSX, "injectSqrtSX", injectSqrtSX );
DebugPrintf ( 1, "done.\n");

%% ----- template search-grid in {f0, tau} ----------
f0Grid  = [min(prior_f0Range):  step_f0 : max(prior_f0Range)];
tauGrid = [min(prior_tauRange): step_tau : max(prior_tauRange)];

%% ----- run 'numSearches' searches
resV = struct();
numSearches = length ( t0V );
for m = 1 : numSearches

  if ( !isempty ( injectionSources ) )
    %% ----- inject QNM signal(s) into analysis segment (wouldn't have changed PSD estimate, which excluded TimeRange)
    multiTS0Inj = injectQNMs ( multiTS0, injectionSources(m) );
  else
    multiTS0Inj = multiTS0;
  endif

  %% ----- narrow-band and whiten/overwhiten data ----------
  [multiTS, multiPSD] = whitenNarrowBandTS ( multiTS0Inj, multiPSD0, FreqRange );

  %% ----- run QNM search(es) over vector of start-times 't0V'
  res_m = searchRingdown ( "multiTS", multiTS, ...
                             "multiPSD", multiPSD, ...
                             "t0", t0V(m), ...
                             "f0Grid", f0Grid, ...
                             "tauGrid", tauGrid, ...
                             "prior_H", prior_H, ...
                             "skyCorr", skyCorr
                           );
  %% 'augment' with derived results for easier plotting
  res_m.f0_est   = credibleInterval ( res_m.posterior_f0, confidence );
  res_m.tau_est  = credibleInterval ( res_m.posterior_tau, confidence );
  res_m.isoConf2 = credibleContourLevel ( res_m.BSG_f0_tau, confidence );

  res_m.tMerger  = tMerger;
  res_m.tOffs    = t0V(m) - tMerger;
  res_m.bname    = sprintf ( "Ringdown-search%04d-GPS%.6fs-f%.0fHz-%.0fHz-tau%.1fms-%.1fms",
                               m, t0V(m), min(prior_f0Range), max(prior_f0Range),
                               1e3 * min(prior_tauRange), 1e3 * max(prior_tauRange)
                             );

  resV(m) = res_m;
  %% ----- save posterior in matrix format ----------
  fname = sprintf ( "%s/%s-BSG.dat", resDir, resV(m).bname );
  tmp = resV(m).BSG_f0_tau;
  save ( "-ascii", fname, "tmp" );

  DebugPrintSummary ( 1, resV(m) );
  %% ---------- plot results summary page
  if ( doPlotSnapshots )
    plotSnapshot ( resV(m), plotMarkers );
    fname = sprintf ( "%s/%s.pdf", resDir, resV(m).bname );
    ezprint ( fname, "width", 512 );
  endif

  %% ----- if injection: quantify PE recovery quality ----------
  if ( !isempty ( injectionSources ) )
    PErecovery(m) = testPErecovery ( resV(m), injectionSources(m) );
  endif %% if injectionSources

endfor %% m = 1 : numSearches

%% ----- save axis {f0, tau} in matrix format ----------
ff0 = resV(1).ff0;
ttau = resV(1).ttau;
save ( "-ascii", strcat(resDir, "/f0s.dat"), "ff0" );
save ( "-ascii", strcat(resDir, "/taus.dat"), "ttau" );

%% ----- plot quantities vs tOffs ----------
if ( doPlotSummary )
  plotSummary ( resV );
  fname = sprintf ( "%s/%s-summary.pdf", resDir, bname );
  ezprint ( fname, "width", 512 );
endif

if ( doPlotContours )
  plotContours ( resV, [], plotMarkers )
  fname = sprintf ( "%s/%s-contours.pdf", resDir, bname );
  ezprint ( fname, "width", 512 );
endif

if ( doPlotH )
  figure(); clf; hold on;
  plot ( [prior_H.x], [prior_H.px], "-xg;prior;", "linestyle", "--", "linewidth", 2 );
  for m = 1 : numSearches
    plot ( [resV(m).posterior_H.x], [resV(m).posterior_H.px], "-b" );
    leg = sprintf ( "+%.1fms", resV(m).tOffs * 1e3 );
    text ( [resV(m).posterior_H.x](end), [resV(m).posterior_H.px](end), leg );
  endfor
  xlabel ( "H" ); ylabel ("pdf");
  grid on; hold off;
  fname = sprintf ( "%s/%s-posterior_H.pdf", resDir, bname );
  ezprint ( fname, "width", 512 );
endif

if ( doPlotBSGHist )
  figure (); clf;
  hist ( [resV.BSG], 100 );
  xlabel ( "<BSG>" );
  fname = sprintf ( "%s/%s-hist.pdf", resDir, bname );
  ezprint ( fname, "width", 512 );
endif

if ( doPlotPErecovery )
  [H_err, H_coverage] = plotPErecovery ( PErecovery );
  set ( 0, "currentfigure", H_err )
  fname = sprintf ( "%s/%s-PE-errors.pdf", resDir, bname );
  ezprint ( fname, "width", 512 );

  set ( 0, "currentfigure", H_coverage )
  fname = sprintf ( "%s/%s-PE-coverage.pdf", resDir, bname );
  ezprint ( fname, "width", 512 );
endif %% doPlotPErecovery()

%% ----- store complete results dump ----------
versioning.octave = version();
versioning.octapps = octapps_gitID();
versioning.scripts = octapps_gitID(".", "QNM-scripts");

fname = sprintf ( "%s/RingdownDriver-%s.hd5", resDir, bname );
save ("-hdf5", fname )
