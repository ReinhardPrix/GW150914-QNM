#!/usr/bin/octave -q

global debugLevel = 1;
%% --------------------------------------------------
%% 1) run a ringdown search as function of start-time 't0' on GW150914
%% --------------------------------------------------

%% ---------- Prior ranges --------------------
tCenter = 1126259462;
tOffsStart = 0.422;		%% this is about when the signal enters f0>~200Hz
dtOffs     = 0.0005;	%% 0.5ms stepsize
tOffsEnd   = 0.438; %% 0.438;

data_FreqRange  = [ 100, 300 ]; %% avoid nasty noise stuff > 300Hz in L1
prior_FreqRange = [ 210, 270 ];
prior_tauRange  = [ 0.2e-3, 20e-3 ];
prior_H         = 4e-22;	%% allow going up from 1e-22 to ~1e-21, fairly "flat" in that region

%% ----- plot data preparation -----
[ts, tsW, tsOW, psd, IFO] = extractTS ( "fMin", min(data_FreqRange), "fMax", max(data_FreqRange), "tCenter", tCenter, "plotResults", false );

%% ----- run search
tOffs = tOffsStart : dtOffs : tOffsEnd;
Nsteps = length(tOffs);

for i = 1:Nsteps
  DebugPrintf ( 1, "tOffs = %.4f s:\n", tOffs(i) );
  ret{i} = searchRingdown ( "tOffs", tOffs(i), "tCenter", tCenter, "prior_FreqRange", prior_FreqRange, "prior_tauRange", prior_tauRange, "prior_H", prior_H, ...
                            "data_FreqRange", data_FreqRange, "plotResults", true );

  %% collect results for plotting
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
  fname = sprintf ( "Results/%s.png", ret{i}.bname );
  ezprint ( fname, "width", 1024, "height", 786, "dpi", 72 );
  %% create symlinks for movie-making
  fnameM = sprintf ( "Results/movie-frame-%02d.png", i - 1 );
  unlink ( fnameM); symlink ( fname, fnameM );
endfor

%% ----- store results dump ----------
fname = sprintf ( "Results/RingdownDriver-%s.hd5", ret{1}.bname );
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

fname = sprintf ( "Results/%s-summary.pdf", ret{1}.bname );
ezprint ( fname, "width", 512 );
