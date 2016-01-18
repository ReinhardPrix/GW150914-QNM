#!/usr/bin/octave -q

global debugLevel = 1;

%% --------------------------------------------------
%% run a ringdown search as function of start-time 't0' on GW150914
%% --------------------------------------------------

%% ========== driver parameters ==========
SFTs = {"./Data/H-1_H1_1800SFT_ER8-C01-1126257832-1800.sft"; "./Data/L-1_L1_1800SFT_ER8-C01-1126258841-1800.sft" };
fSamp = 4084;	%% chosen such that 7.1ms timeshift can be represented by ~integer bin-shift: 7.1e-3 * 4084 =  28.9964000000000

plotSummary = false;
plotSpectra = false;
useTSBuffer = true;
plotPosteriors = true;

if ( !exist ("searchType" ) )
  searchType = "verify";
endif
if ( !exist ( "extraLabel" ) )
  extraLabel = "";
endif

%% ---------- Prior range defaults ----------
data_FreqRange  = [ 100, 300 ]; %% avoid nasty noise stuff > 300Hz in L1
prior_f0Range = [ 210, 270 ];
prior_tauRange  = [ 1e-3, 20e-3 ];
prior_H         = 4e-22;	%% allow going up from 1e-22 to ~1e-21, fairly "flat" in that region
step_f0 = 0.1;
step_tau = 0.2e-3;

switch ( searchType )
  case "verify"
    %% ---------- test-case to compare different code-versions on ----------
    step_f0 = 0.5;
    step_tau = 0.5e-3;
    tCenter = 1126259462;
    if ( !exist ( "tOffs" ) )
      tOffs = 0.43;
    endif
    tOffsStart = tOffs;
    dtOffs     = 0.0005;
    tOffsEnd   = tOffs;
    sideband = 15;
    data_FreqRange  = [ 15 + sideband, 2000 - sideband - 1 ];
    plotSpectra = true;
    useTSBuffer = true;
    plotBSGHist = false;

  case "onSource"
    %% ---------- "ON-SOURCE ----------
    tCenter = 1126259462;
    tOffsStart = 0.422;		%% this is about when the signal enters f0>~200Hz
    dtOffs     = 0.0005;	%% 0.5ms stepsize
    tOffsEnd   = 0.438; %% 0.438;
    plotBSGHist = false;
    plotSummary = true;
    plotPosteriors = true;

  case "offSource"
    %% ---------- "OFF-SOURCE" for background estimation ----------
    tCenter     = 1126259472; %% 10s past GW150914
    tOffsStart  = -5;
    dtOffs      = 0.1;	%% 100ms stepsize, to avoid template overlap --> 'independent' templates
    tOffsEnd    = 5;
    plotPosteriors = true;
    plotSummary = true;
    plotSpectra = false;
    useTSBuffer = true;
    plotBSGHist = true;
  otherwise
    error ("Unknown searchType = '%s' specified\n", searchType );

endswitch

%% ----- data preparation -----

%% load frequency-domain data from SFTs:
for X = 1:length(SFTs)
  [ts{X}, psd{X}] = extractTSfromSFT ( "SFTpath", SFTs{X}, "fMin", min(data_FreqRange), "fMax", max(data_FreqRange), "fSamp", fSamp, "tCenter", tCenter, "plotSpectrum", plotSpectra, "useBuffer", useTSBuffer );
endfor

%% create unique time-tagged 'ResultsDir' for each run:
gm = gmtime ( time () );
resDir = sprintf ( "Results/Results-%02d%02d%02d-%02dh%02d-%s%s", gm.year - 100, gm.mon + 1, gm.mday, gm.hour, gm.min, searchType, extraLabel );
[status, msg, id] = mkdir ( resDir ); assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );
addpath ( pwd() );
cd ( resDir );

try
%% ----- run search
tOffs = tOffsStart : dtOffs : tOffsEnd;
Nsteps = length(tOffs);

ret = cell ( 1, Nsteps );
for i = 1:Nsteps
  DebugPrintf ( 1, "tOffs = %.4f s:\n", tOffs(i) );
  ret{i} = searchRingdown ( "ts", ts, "psd", psd, "tOffs", tOffs(i), "tCenter", tCenter, ...
                            "prior_f0Range", prior_f0Range, "step_f0", step_f0, ...
                            "prior_tauRange", prior_tauRange, "step_tau", step_tau, ...
                            "prior_H", prior_H, ...
                            "plotResults", plotPosteriors
                          );

  %% ----- save posterior in matrix format ----------
  if ( i == 1 )
    fname = sprintf ( "f0s.dat", ret{i}.bname );
    tmp = ret{i}.ff0;
    save ( "-ascii", fname, "tmp" );
    fname = sprintf ( "taus.dat", ret{i}.bname );
    tmp = ret{i}.ttau;
    save ( "-ascii", fname, "tmp" );
  endif
  fname = sprintf ( "%s-BSG.dat", ret{i}.bname );
  tmp = ret{i}.BSG;
  save ( "-ascii", fname, "tmp" );

  fname = sprintf ( "%s-SNR.dat", ret{i}.bname );
  tmp = ret{i}.SNR;
  save ( "-ascii", fname, "tmp" );

  %% ---------- collect results for plotting
  BSG_mean(i) = ret{i}.BSG_mean;
  A_MPE(i)    = ret{i}.A_MPE;
  phi0_MPE(i) = ret{i}.phi0_MPE;
  SNR_MPE(i)  = ret{i}.SNR_MPE;

  f0_MPE(i)   = ret{i}.f0_est.MPE;
  f0_lerr(i)  = f0_MPE(i) - ret{i}.f0_est.lower;
  f0_uerr(i)  = ret{i}.f0_est.upper - f0_MPE(i);

  tau_MPE(i)   = ret{i}.tau_est.MPE;
  tau_lerr(i)  = tau_MPE(i) - ret{i}.tau_est.lower;
  tau_uerr(i)  = ret{i}.tau_est.upper - tau_MPE(i);

  %% ---------- plot results summary page
  %%fname = sprintf ( "%s.png", ret{i}.bname );
  %%ezprint ( fname, "width", 1024, "height", 786, "dpi", 72 );
  if ( plotPosteriors )
    fname = sprintf ( "%s.pdf", ret{i}.bname );
    ezprint ( fname, "width", 512 );
  endif

endfor

%% ----- store results dump ----------
fname = sprintf ( "RingdownDriver-%s.hd5", ret{1}.bname );
save ("-hdf5", fname )

%% ----- plot quantities vs tOffs ----------
if ( plotSummary )
  figure(); clf;
  xrange = [ tOffsStart, tOffsEnd ];

  subplot ( 3, 1, 1, "align" );
  semilogy ( tOffs, BSG_mean, "-o" ); grid on;
  xlim ( xrange );
  ylabel ("<BSG>");

  subplot ( 3, 1, 2, "align" );
  errorbar ( tOffs, f0_MPE, f0_lerr, f0_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( prior_f0Range );
  ylabel ("f0 [Hz]");

  subplot ( 3, 1, 3, "align" );
  errorbar ( tOffs, 1e3*tau_MPE, 1e3*tau_lerr, 1e3*tau_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylabel ("tau [ms]");
  xlabel ( sprintf ( "%.0f + tOffs [s]", tCenter) );
  fname = sprintf ( "%s-summary.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );
endif

if ( plotBSGHist )
  figure(); clf;
  hist ( BSG_mean, 20 );
  xlabel ( "<BSG>" );
  fname = sprintf ( "%s-hist.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );
endif

cd ("../..");

catch
  err = lasterror();
  warning(err.identifier, err.message);
  err.stack.file
  err.stack.name
  err.stack.line

  cd ("../..");	%% make sure we end up in main dir in case of failure
end_try_catch
