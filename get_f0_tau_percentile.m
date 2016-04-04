function percentile = get_f0_tau_percentile ( f0Val, tauVal, ff0, ttau, BSG_f0_tau )
  %% determine percentile of given values {f0Val, tauVal} in the 2D posterior {ff0, ttau, BSG_f0_tau}
  %% where percentiles are measured at iso-probability 'heights' from the peak (MPE)

  assert ( isscalar ( f0Val ) && isscalar ( tauVal ) );
  assert ( (size ( ff0 ) == size ( ttau )) && (size(ff0) == size(BSG_f0_tau)) );

  f0s = unique ( ff0(:) );
  taus = unique ( ttau(:) );

  [x, i_f0]  = min ( abs ( f0Val - f0s(:) ) );
  [x, i_tau] = min ( abs ( tauVal - taus(:) ) );

  height_val = BSG_f0_tau ( i_tau, i_f0 );
  total_post = sum ( BSG_f0_tau(:) );
  above_post = sum ( BSG_f0_tau ( BSG_f0_tau >= height_val ) );

  percentile = above_post / total_post;

  return;

endfunction %% get_f0_tau_percentile()
