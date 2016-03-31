function DebugPrintSummary ( level, res_l, resCommon )

  %% summarize numerical outcomes on stdout
  DebugPrintf (level, "---------- t0 = tMerger + tOffs = %.6f s + %.1f ms = %.6f s ----------\n", resCommon.tMerger, res_l.tOffs * 1e3, res_l.t0 );
  DebugPrintf (level, "log10<BSG>= %.2g\n", log10(res_l.BSG) );
  DebugPrintf (level, "f0_est    = {%.1f -%.1f +%.1f} Hz\n", res_l.f0_est.MPE, res_l.f0_est.lerr, res_l.f0_est.uerr );
  DebugPrintf (level, "tau_est   = {%.1f -%.1f +%1.f} ms\n", 1e3 * [res_l.tau_est.MPE, res_l.tau_est.lerr, res_l.tau_est.uerr] );
  DebugPrintf (level, "-\n");
  DebugPrintf (level, "A_MP      = %.2g\n", res_l.AmpMP.A );
  DebugPrintf (level, "phi0_MP   = %.2g\n", res_l.AmpMP.phi0 );
  DebugPrintf (level, "SNR_MP    = %.2g\n", res_l.AmpMP.SNR );
  DebugPrintf (level, "-\n");
  DebugPrintf (level, "A_ML      = %.2g\n", res_l.AmpML.A );
  DebugPrintf (level, "phi0_ML   = %.2g\n", res_l.AmpML.phi0 );
  DebugPrintf (level, "SNR_ML    = %.2g\n", res_l.AmpML.SNR );
  DebugPrintf (level, "-\n");
  return;
endfunction
