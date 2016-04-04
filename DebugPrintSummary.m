function DebugPrintSummary ( level, res )

  %% summarize numerical outcomes on stdout
  DebugPrintf (level, "---------- t0 = tMerger + tOffs = %.6f s + %.1f ms = %.6f s ----------\n", res.tMerger, res.tOffs * 1e3, res.t0 );
  DebugPrintf (level, "log10<BSG>= %.2g\n", log10(res.BSG) );
  if ( !isempty ( res.f0_est ) )
    DebugPrintf (level, "f0_est    = {%.1f -%.1f +%.1f} Hz\n", res.f0_est.MPE, res.f0_est.lerr, res.f0_est.uerr );
  endif
  if ( !isempty ( res.tau_est ) )
    DebugPrintf (level, "tau_est   = {%.1f -%.1f +%1.f} ms\n", 1e3 * [res.tau_est.MPE, res.tau_est.lerr, res.tau_est.uerr] );
  endif
  DebugPrintf (level, "A    [MP|ML] = %.2g | %.2g\n", res.AmpMP.A, res.AmpML.A );
  DebugPrintf (level, "phi0 [MP|ML] = %.2g | %.2g\n", res.AmpMP.phi0, res.AmpML.phi0 );
  DebugPrintf (level, "SNR  [MP|ML] = %.2g | %.2g\n", res.AmpMP.SNR, res.AmpML.SNR );
  return;
endfunction
