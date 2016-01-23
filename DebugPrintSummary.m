function DebugPrintSummary ( level, in )
  global tMergerOffs;
  global tEvent;

  %% summarize numerical outcomes on stdout
  DebugPrintf (level, "tGPS = %.0f + %f s = tMerger + %.2fms\n", fix(in.tGPS), rem(in.tGPS,1), (in.tGPS - (tEvent + tMergerOffs)) * 1e3 );
  DebugPrintf (level, "log10<BSG>    = %.2g\n", log10(in.BSG_mean) );
  DebugPrintf (level, "f0_est  = {%.1f -%.1f +%.1f} Hz\n", in.f0_est.MPE, in.f0_est.lerr, in.f0_est.uerr );
  DebugPrintf (level, "tau_est = {%.1f -%.1f +%1.f} ms\n", 1e3 * [in.tau_est.MPE, in.tau_est.lerr, in.tau_est.uerr] );
  DebugPrintf (level, "A_MPE   = %.2g\n", in.A_MPE2 );
  DebugPrintf (level, "phi0_MPE= %.2g\n", in.phi0_MPE2 );
  DebugPrintf (level, "SNR_MPE = %.2g\n", in.SNR_MPE2 );

endfunction
