function PErecovery = testPErecovery ( res, injectionSource )
  %% quantify PE recovery compared to actual injection

  ff0  = res.ff0;
  ttau = res.ttau;
  f0   = unique ( ff0 );
  tau  = unique ( ttau );

  PErecovery = struct();

  A_inj     = injectionSource.A;
  phi0_inj  = injectionSource.phi0;
  f0_inj    = injectionSource.f0;
  tau_inj   = injectionSource.tau;

  AmpEst    = res.AmpMP;

  PErecovery.t0_offs  = (res.t0 - injectionSource.t0);

  PErecovery.A_relerr   = (AmpEst.A - A_inj) / A_inj;

  Dphi0   = (AmpEst.phi0 - phi0_inj);
  PErecovery.phi0_relerr = min ( abs ( rem ( [Dphi0, Dphi0 + 2*pi], 2*pi ) ) ) / (2*pi);
  PErecovery.f0_relerr   = (res.lambdaMP.f0 - f0_inj) / f0_inj;
  PErecovery.tau_relerr  = (res.lambdaMP.tau - tau_inj) / tau_inj;

  A_s = -A_inj * sin(phi0_inj);
  A_c =  A_inj * cos(phi0_inj);
  [x, i_f0]  = min ( abs ( f0_inj - f0(:) ) );
  [x, i_tau] = min ( abs ( tau_inj - tau(:) ) );
  M_ss = res.Mxy.ss( i_tau, i_f0 );
  M_sc = res.Mxy.sc( i_tau, i_f0 );
  M_cc = res.Mxy.cc( i_tau, i_f0 );
  SNR_inj = sqrt ( A_s^2 * M_ss + 2 * A_s * M_sc * A_c + A_c^2 * M_cc );

  PErecovery.SNR_inj   = SNR_inj;
  PErecovery.SNR_ML    = res.AmpML.SNR;
  PErecovery.SNR_MP    = res.AmpMP.SNR;
  PErecovery.f0_tau_percentile = get_f0_tau_percentile ( f0_inj, tau_inj, ff0, ttau, res.BSG_f0_tau );
  PErecovery.BSG       = res.BSG;

  return;

endfunction %% testPErecovery()
