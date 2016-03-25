function DebugPrintSummary ( level, res_l, resCommon )

  %% summarize numerical outcomes on stdout
  DebugPrintf (level, "t0 = tMerger + tOffs = %.6f s + %.1f ms = %.6f s\n", resCommon.tMerger, res_l.tOffs * 1e3, res_l.t0 );
  DebugPrintf (level, "log10<BSG>= %.2g\n", log10(res_l.BSG) );
  DebugPrintf (level, "f0_est    = {%.1f -%.1f +%.1f} Hz\n", res_l.f0_est.MPE, res_l.f0_est.lerr, res_l.f0_est.uerr );
  DebugPrintf (level, "tau_est   = {%.1f -%.1f +%1.f} ms\n", 1e3 * [res_l.tau_est.MPE, res_l.tau_est.lerr, res_l.tau_est.uerr] );
  DebugPrintf (level, "A_MP2D    = %.2g\n", res_l.A_MP2D );
  DebugPrintf (level, "phi0_MP2D = %.2g\n", res_l.phi0_MP2D );
  DebugPrintf (level, "SNR_MP2D  = %.2g\n", res_l.SNR_MP2D );
  DebugPrintf (level, "H_MPE     = %.2g\n", res_l.H_MPE );

endfunction
