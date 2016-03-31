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
shiftL1 = 7.0e-3;	%% time-shift to apply to L1 data-stream: currently 'official' value (v8)

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
if ( !exist("psd_version") )    global psd_version = 2; endif
if ( !exist("cleanLines") )     global cleanLines = false; endif
if ( !exist("data_FreqRange") ) data_FreqRange  = [ 30, 1e3 ]; endif
if ( !exist("injectionSources") ) injectionSources = []; endif
if ( !exist("assumeSqrtSX") ) 	assumeSqrtSX = cell(1,numIFOs); endif	%% empty by default
if ( !exist("injectSqrtSX") ) 	injectSqrtSX = cell(1,numIFOs); endif	%% empty by default
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

    doPlotContours  = true;
    doPlotSummary   = false;
    doPlotSnapshots = true;
    doPlotSpectra   = true;
    doPlotH         = true;

  case "onSource"
    %% ---------- "ON-SOURCE ----------
    tMerger = tMergerGW150914;
    tOffsV = [ 1, 3, 5, 7 ] * 1e-3;

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

    doPlotSnapshots = false;
    doPlotContours  = false;
    doPlotSummary   = true;
    doPlotSpectra   = false;
    doPlotBSGHist   = true;
    doPlotH         = false;

  case "injections"
    %% NOTE: injections defaults to injections into real-data, but override with 'injectSqrtSX' and 'assumeSqrtSX' settings
    tMerger = tMergerGW150914 + 10;;	%% go to some 'off source' time-stretch
    tOffsV = [ -3.0 : 0.2 : 3.0 ];
    clear("injectionSources");
    for l = 1 : length(tOffsV)
      if ( noMismatchInj )
        Nf0  = floor ( diff(prior_f0Range) / step_f0 ) + 1;
        Ntau = floor ( diff(prior_tauRange) / step_tau ) + 1;
        i_f0 = unidrnd ( Nf0 );
        i_tau = unidrnd ( Ntau );
        f0_inj  = min(prior_f0Range) + (i_f0 - 1) * step_f0;
        tau_inj = min(prior_tauRange) + (i_tau - 1) * step_tau;
      else
        f0_inj = unifrnd ( prior_f0Range(1), prior_f0Range(2) );
        tau_inj = unifrnd ( prior_tauRange(1), prior_tauRange(2) );
      endif
      injectionSources(l) = struct ( "name", 	sprintf("Inj-%d", l), ...
                                     "t0", 	tMerger + tOffsV(l), ...
                                     "A", 	unifrnd ( 3e-22, 8e-22 ), ...
                                     "phi0", 	unifrnd ( 0, 2*pi ), ...
                                     "f0", 	f0_inj, ...
                                     "tau", 	tau_inj, ...
                                     "shiftL1", shiftL1 ...
                                   );
    endfor %% i = l:numInjections

    data_FreqRange  = [ 30, 1e3 ];	%% full data-range, better PE-recovery at low SNR(?)
    doPlotContours  = false;
    doPlotSummary   = false;
    doPlotSnapshots = false;
    doPlotSpectra   = false;
    doPlotH         = false;
    doPlotPErecovery= true;

  otherwise
    error ("Unknown searchType = '%s' specified\n", searchType );

endswitch

%% ----- data preparation -----

%% load frequency-domain data from SFTs:
for X = 1:length(SFTs)
  [ts{X}, psd{X}] = extractTSfromSFT ( "SFTpath", SFTs{X}, ...
                                       "fMin", min(data_FreqRange), ...
                                       "fMax", max(data_FreqRange), ...
                                       "fSamp", fSamp, ...
                                       "tCenter", fix(tMerger), ...
                                       "Twindow", 4, ...
                                       "plotSpectrum", doPlotSpectra, ...
                                       "injectionSources", injectionSources, ...
                                       "assumeSqrtSn", assumeSqrtSX{X}, ...
                                       "injectSqrtSn", injectSqrtSX{X} ...
                                     );
  %% for plotting OW-timeseries: store noise-values at 'GR frequency f0GR'
  [val, freqInd] = min ( abs ( psd{X}.fk - f0GR.val ) );
  ts{X}.SX_GR = psd{X}.Sn ( freqInd );
  DebugPrintf ( 1, "X = %d: sqrt(SX) = %g /sqrt(Hz)\n", X, sqrt ( ts{X}.SX_GR ) );
endfor




%% create unique time-tagged 'ResultsDir' for each run:
gm = gmtime ( time () );
resDir = sprintf ( "Results/Results-%02d%02d%02d-%02dh%02d-%s-data%.0fHz-%.0fHz-psd_v%d-lineCleaning%s%s",
                   gm.year - 100, gm.mon + 1, gm.mday, gm.hour, gm.min, searchType, data_FreqRange,
                   psd_version, ifelse ( cleanLines, "On", "Off" ),
                   extraLabel );
if ( noMismatchInj )
  resDir = strcat ( resDir, "-noMismatchInj" );
endif
if ( !isempty ( injectSqrtSX{1} ) )
  resDir = strcat ( resDir, "-injectSqrtSX", sprintf("_%.2g", cell2mat(injectSqrtSX) ) );
endif
if ( !isempty ( assumeSqrtSX{1} ) )
  resDir = strcat ( resDir, "-assumeSqrtSX", sprintf("_%.2g", cell2mat(assumeSqrtSX) ) );
endif

[status, msg, id] = mkdir ( resDir ); assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );
addpath ( pwd() );
cd ( resDir );

%%try
%% ----- run search
Nsearches = length(tOffsV);
PErecovery = [];
clear ("res" );
[resV, resCommon] = searchRingdown ( "ts", ts, "psd", psd, ...
                                     "t0V", tMerger + tOffsV, ...
                                     "prior_f0Range", prior_f0Range, "step_f0", step_f0, ...
                                     "prior_tauRange", prior_tauRange, "step_tau", step_tau, ...
                                     "prior_H", prior_H, ...
                                     "shiftL1", shiftL1
                                   );
assert ( length(resV) == Nsearches );
resCommon.tMerger = tMerger;

for l = 1 : Nsearches
  resV(l).tOffs = resV(l).t0 - resCommon.tMerger;

  %% ----- save posterior in matrix format ----------
  fname = sprintf ( "%s-BSG.dat", resV(l).bname );
  tmp = resV(l).BSG_f0_tau;
  save ( "-ascii", fname, "tmp" );

  %% ----- compute derived quantities
  %% 1D marginalized posteriors on {f,tau} ----------
  posterior2D = resV(l).BSG_f0_tau;
  tmp = sum ( posterior2D, 1 );
  resV(l).posterior_f0   = tmp / sum ( tmp(:) );
  tmp = sum ( posterior2D, 2 );
  resV(l).posterior_tau   = tmp / sum ( tmp(:) );

  ff0 = resCommon.ff0;
  ttau = resCommon.ttau;
  f0  = unique ( ff0 );
  tau = unique ( ttau );
  resV(l).f0_est   = credibleInterval ( f0,   resV(l).posterior_f0,   confidence );
  resV(l).tau_est  = credibleInterval ( tau,  resV(l).posterior_tau,  confidence );
  resV(l).isoConf2 = credibleContourLevel ( posterior2D, confidence );

  DebugPrintSummary ( 1, resV(l), resCommon );
  %% ---------- plot results summary page
  if ( doPlotSnapshots )
    plotSnapshot ( resV(l), resCommon, plotMarkers );
  endif

  %% ----- if injection: quantify PE recovery quality ----------
  if ( !isempty ( injectionSources ) )

    %% single-QNM search: compare to l'th injection only, assuming we're targeting one injection per 'step'
    A_inj = injectionSources(l).A;
    phi0_inj = injectionSources(l).phi0;
    f0_inj = injectionSources(l).f0;
    tau_inj = injectionSources(l).tau;
    AmpEst = resV(l).AmpMP;

    PErecovery(end+1).t0_offs  = (resV(l).t0 - injectionSources(l).t0);

    PErecovery(end).A_relerr   = (AmpEst.A - A_inj) / A_inj;
    Dphi0   = (AmpEst.phi0 - phi0_inj);
    PErecovery(end).phi0_relerr = min ( abs ( rem ( [Dphi0, Dphi0 + 2*pi], 2*pi ) ) ) / (2*pi);
    PErecovery(end).f0_relerr  = (resV(l).lambdaMP.f0 - f0_inj) / f0_inj;
    PErecovery(end).tau_relerr = (resV(l).lambdaMP.tau - tau_inj) / tau_inj;

    A_s = -A_inj * sin(phi0_inj);
    A_c =  A_inj * cos(phi0_inj);
    [x, i_f0]  = min ( abs ( f0_inj - f0(:) ) );
    [x, i_tau] = min ( abs ( tau_inj - tau(:) ) );
    M_ss = resCommon.Mxy.ss( i_tau, i_f0 );
    M_sc = resCommon.Mxy.sc( i_tau, i_f0 );
    M_cc = resCommon.Mxy.cc( i_tau, i_f0 );
    SNR_inj = sqrt ( A_s^2 * M_ss + 2 * A_s * M_sc * A_c + A_c^2 * M_cc );

    PErecovery(end).SNR_inj   = SNR_inj;
    PErecovery(end).SNR_ML    = resV(l).AmpML.SNR;
    PErecovery(end).SNR_MP    = resV(l).AmpMP.SNR;
    PErecovery(end).f0_tau_percentile = get_f0_tau_percentile ( f0_inj, tau_inj, ff0, ttau, posterior2D );
    PErecovery(end).BSG       = resV(l).BSG;
    plotSnapshot ( resV(l), resCommon, injectionSources(l) );
  endif %% if injectionSources

endfor %% l = 1:Nsearches

%% ----- save axis {f0, tau} in matrix format ----------
save ( "-ascii", "f0s.dat", "ff0" );
save ( "-ascii", "taus.dat", "ttau" );

%% ----- store complete results dump ----------
versioning.octave = version();
versioning.octapps = octapps_gitID();
versioning.scripts = octapps_gitID(".", "QNM-scripts");

fname = sprintf ( "RingdownDriver-%s.hd5", resCommon.bname );
save ("-hdf5", fname )

%% ----- plot quantities vs tOffs ----------
if ( doPlotSummary )
  plotSummary ( resV, resCommon );
endif

if ( doPlotContours )
  plotContours ( resV, resCommon, [], plotMarkers )
endif

if ( doPlotH )
  figure(); clf; hold on;
  plot ( prior_H(:,1), prior_H(:,2), "-xg;prior;", "linestyle", "--", "linewidth", 2 );
  for i = 1 : Nsearches
    plot ( prior_H(:,1), resV(l).post_H, "-b" );
    leg = sprintf ( "+%.1fms", resV(l).tOffs * 1e3 );
    text ( prior_H(end,1), resV(l).post_H(end), leg );
  endfor
  xlabel ( "H" ); ylabel ("pdf");
  grid on; hold off;
  fname = sprintf ( "%s-post_H.pdf", resCommon.bname );
  ezprint ( fname, "width", 512 );
endif

if ( doPlotBSGHist )
  figure (); clf;
  hist ( [resV.BSG], 100 );
  xlabel ( "<BSG>" );
  fname = sprintf ( "%s-hist.pdf", resCommon.bname );
  ezprint ( fname, "width", 512 );
endif

if ( doPlotPErecovery )
  plotPErecovery ( PErecovery, resCommon.bname );
endif %% doPlotPErecovery()

cd ("../..");

## catch
##   err = lasterror();
##   warning(err.identifier, err.message);
##   err.stack.file
##   err.stack.name
##   err.stack.line

##   cd ("../..");	%% make sure we end up in main dir in case of failure
## end_try_catch


