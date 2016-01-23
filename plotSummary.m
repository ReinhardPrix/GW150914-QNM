function plotSummary ( in )
  global iFig0 = 0;
  global tMergerOffs
  global tEvent;
  global f0GR;
  global taumsGR;

  %% ----- prepare quantities to be plotted versus QNM start-time
  Nsteps = length ( in );
  tOffs = BSG_mean = SNR_MPE2 = f0_MPE = f0_lerr = f0_uerr = taums_MPE = taums_lerr = taums_uerr = zeros ( 1, Nsteps );
  for i = 1 : Nsteps
    tOffs(i)      = in{i}.tGPS - tEvent;
    BSG_mean(i)   = in{i}.BSG_mean;
    SNR_MPE2(i)   = in{i}.SNR_MPE2;

    f0_MPE(i)     = in{i}.f0_est.MPE;
    f0_lerr(i)    = in{i}.f0_est.lerr;
    f0_uerr(i)    = in{i}.f0_est.uerr;

    taums_MPE(i)  = in{i}.tau_est.MPE * 1e3;
    taums_lerr(i) = in{i}.tau_est.lerr * 1e3;
    taums_uerr(i) = in{i}.tau_est.uerr * 1e3;
  endfor

  tOffs_Range = [ (min ( tOffs(:) )), (max ( tOffs(:) )) ];
  taums_Range = [ (min ( in{1}.ttau(:) )), (max ( in{1}.ttau(:))) ] * 1e3;
  f0_Range    = [ (min ( in{1}.ff0(:) )),  (max ( in{1}.ff0(:) )) ];

  %% ===== plot summary 1 ==========
  figure ( iFig0 + 4 ); clf;

  %% ----- plot log10BSG(tOffs)
  subplot ( 2, 2, 1, "align" );
  xrange = tOffs_Range;
  yrange = [ -1, 10 ];
  plot ( tOffs, log10(BSG_mean), "-o" );
  xlim ( xrange );
  ylim ( yrange );
  line ( xrange, 0, "linestyle", "-", "linewidth", 3 );
  line ( xrange, 1, "linestyle", ":", "linewidth", 3 );
  line ( tMergerOffs * [1,1], yrange, "linestyle", "-", "linewidth", 2 );
  xlim ( xrange );
  ylabel ("log10<BSG>");
  xlabel ( "tOffs [s]")
  grid on;
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );
  xlabel ( "tOffs from Merger [ms]")

  %% ----- plot f0_MPE(tOffs)
  subplot ( 2, 2, 2, "align" );
  errorbar ( tOffs, f0_MPE, f0_lerr, f0_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( f0_Range );
  ylabel ("f0 [Hz]");
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );
  xlabel ( "tOffs from Merger [ms]")


  %% ----- plot SNR(tOffs)
  subplot ( 2, 2, 3, "align" );
  plot ( tOffs, SNR_MPE2, "-o" ); grid on;
  xlim ( xrange );
  ylabel ("SNR(MPE)");
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );


  %% ----- plot tau_MPE(tOffs)
  subplot ( 2, 2, 4, "align" );
  errorbar ( tOffs, taums_MPE, taums_lerr, taums_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( taums_Range );
  ylabel ("tau [ms]");
  xlabel ( sprintf ( "%.0f + tOffs [s]", tEvent) );
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );

  fname = sprintf ( "%s-summary.pdf", in{1}.bname );
  ezprint ( fname, "width", 512 );


  return;

endfunction
