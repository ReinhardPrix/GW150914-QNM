function DebugPrintSummary ( level, res )

  %% summarize numerical outcomes on stdout
  DebugPrintf (level, "---------- t0 = tMerger + tOffs = %.6f s + %.1f ms = %.6f s ----------\n", res.tMerger, res.tOffs * 1e3, res.t0 );
  DebugPrintf (level, "log10<BSG>= %.2g\n", log10(res.BSG) );
  DebugPrintf (level, "f0_est    = {%.1f -%.1f +%.1f} Hz\n", res.f0_est.MPE, res.f0_est.lerr, res.f0_est.uerr );
  DebugPrintf (level, "tau_est   = {%.1f -%.1f +%1.f} ms\n", 1e3 * [res.tau_est.MPE, res.tau_est.lerr, res.tau_est.uerr] );
  DebugPrintf (level, "-\n");
  DebugPrintf (level, "A_MP      = %.2g\n", res.AmpMP.A );
  DebugPrintf (level, "phi0_MP   = %.2g\n", res.AmpMP.phi0 );
  DebugPrintf (level, "SNR_MP    = %.2g\n", res.AmpMP.SNR );
  DebugPrintf (level, "-\n");
  DebugPrintf (level, "A_ML      = %.2g\n", res.AmpML.A );
  DebugPrintf (level, "phi0_ML   = %.2g\n", res.AmpML.phi0 );
  DebugPrintf (level, "SNR_ML    = %.2g\n", res.AmpML.SNR );
  DebugPrintf (level, "-\n");
  return;
endfunction
