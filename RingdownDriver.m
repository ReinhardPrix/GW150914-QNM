#!/usr/bin/octave -q

global debugLevel = 1;

%% --------------------------------------------------
%% run a ringdown search as function of start-time 't0' on GW150914
%% --------------------------------------------------

%% ========== driver parameters ==========
onSource = false;

%% ========== Prior ranges ==========

%% ---------- "ON-SOURCE RANGE" ----------
if ( onSource )
  tCenter = 1126259462;
  tOffsStart = 0.422;		%% this is about when the signal enters f0>~200Hz
  dtOffs     = 0.0005;	%% 0.5ms stepsize
  tOffsEnd   = 0.438; %% 0.438;
  extraLabel = "";	%% provide extra info about what's specific in this run
else
  %% ---------- "OFF-SOURCE RANGE" for background estimation ----------
  tCenter = 1126259472;
  tOffsStart = -1;		%% this is about when the signal enters f0>~200Hz
  dtOffs     = 0.1;	%% 100ms stepsize, to avoid template overlap --> 'independent' templates
  tOffsEnd   = 1;
  extraLabel = "-OffSource";
endif

data_FreqRange  = [ 100, 300 ]; %% avoid nasty noise stuff > 300Hz in L1
prior_FreqRange = [ 210, 270 ];
prior_tauRange  = [ 0.2e-3, 20e-3 ];
prior_H         = 4e-22;	%% allow going up from 1e-22 to ~1e-21, fairly "flat" in that region

%% ----- data preparation -----
DebugPrintf ( 1, "Extracting timeseries ... ");
[ts, tsW, tsOW, psd] = extractTS ( "fMin", min(data_FreqRange), "fMax", max(data_FreqRange), "tCenter", tCenter, "plotResults", false );
DebugPrintf ( 1, "done.\n");

%% create unique time-tagged 'ResultsDir' for each run:
gm = gmtime ( time () );
resDir = sprintf ( "Results-%02d%02d%02d-%02dh%02d%s", gm.year - 100, gm.mon, gm.mday, gm.hour, gm.min, extraLabel );
[status, msg, id] = mkdir ( resDir ); assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );
addpath ( pwd() );
cd ( resDir );

try
%% ----- run search
tOffs = tOffsStart : dtOffs : tOffsEnd;
Nsteps = length(tOffs);

for i = 1:Nsteps
  DebugPrintf ( 1, "tOffs = %.4f s:\n", tOffs(i) );
  ret{i} = searchRingdown ( "tsOW", tsOW, "psd", psd, "tOffs", tOffs(i), "tCenter", tCenter, "prior_FreqRange", prior_FreqRange, "prior_tauRange", prior_tauRange, "prior_H", prior_H, "plotResults", true );

  %% ----- save posterior in matrix format ----------
  if ( i == 1 )
    fname = sprintf ( "Freqs.dat", ret{i}.bname );
    tmp = ret{i}.posterior.Freq;
    save ( "-ascii", fname, "tmp" );
    fname = sprintf ( "taus.dat", ret{i}.bname );
    tmp = ret{i}.posterior.tau;
    save ( "-ascii", fname, "tmp" );
  endif
  fname = sprintf ( "%s-BSG.dat", ret{i}.bname );
  tmp = ret{i}.posterior.BSG;
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
  fname = sprintf ( "%s.pdf", ret{i}.bname );
  ezprint ( fname, "width", 512 );

endfor

%% ----- store results dump ----------
fname = sprintf ( "RingdownDriver-%s.hd5", ret{1}.bname );
save ("-hdf5", fname )

%% ----- plot quantities vs tOffs ----------
figure(); clf;

subplot ( 3, 1, 1, "align" );
semilogy ( tOffs, BSG_mean, "-o" ); grid on;
yrange = ylim(); ylim ( [ 1e-1, max(yrange) ] );
ylabel ("<BSG>");

subplot ( 3, 1, 2, "align" );
errorbar ( tOffs, f0_MPE, f0_lerr, f0_uerr, ";90%;" ); grid on;
ylim ( prior_FreqRange );
ylabel ("f0 [Hz]");

subplot ( 3, 1, 3, "align" );
errorbar ( tOffs, 1e3*tau_MPE, 1e3*tau_lerr, 1e3*tau_uerr, ";90%;" ); grid on;
ylabel ("tau [ms]");
xlabel ("tOffs [s]");

fname = sprintf ( "%s-summary.pdf", ret{1}.bname );
ezprint ( fname, "width", 512 );

cd ("..");

catch err
  err
  lasterror ( err );
  cd ("..");	%% make sure we end up in main dir in case of failure
end_try_catch
