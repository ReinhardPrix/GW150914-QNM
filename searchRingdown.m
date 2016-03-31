## Copyright (C) 2015, 2016 Reinhard Prix
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

function [ resV, resCommon ] = searchRingdown ( varargin )
  global debugLevel = 1;
  global psd_version = 1;
  global cleanLines = false;

  uvar = parseOptions ( varargin,
                        {"ts", "cell" },	%% cell-array [over detectors]: normal, whitentend, and over-whitened timeseries
                        {"psd", "cell"},	%% cell-array [over detectors]: PSD estimate over frequency range, for each detector
                        {"t0V", "real,strictpos,vector", 1126259462.43 },	%% ringdown start-time in GPS s
                        {"prior_f0Range", "real,strictpos,vector", [220,  270] },
                        {"step_f0", "real,strictpos,scalar", 0.1 },
                        {"prior_tauRange", "real,vector", [1e-3, 30e-3] },
                        {"step_tau", "real,strictpos,scalar", 0.2e-3 },
                        {"prior_H", "real,strictpos,matrix", [4e-22, 1]},
                        {"plotResults", "bool", false },
                        {"shiftL1", "real,scalar", 7.0e-3}	%% time-shift to apply to L1 data-stream: currently 'official' value (v8)
                      );

  if ( isscalar ( uvar.prior_H ) )
    numH = 1;
    priorH = [ uvar.prior_H, 1 ];
    priorHname = sprintf ( "%.2g", uvar.prior_H );
  else
    priorHname = "Jeffreys";
    [numH, check] = size ( uvar.prior_H );
    assert ( check == 2, "prior_H needs to be of size [numH x 2]\n");
    priorH = uvar.prior_H;
    normH = sum ( priorH ( :, 2 ) );
    priorH ( :, 2) /= normH;
  endif

  psd = uvar.psd;
  Ndet = length(uvar.ts);
  assert ( Ndet == 2, "Currently only 2 detectors supported, must be 'H1' and 'L1'\n");
  fk = psd{1}.fk;

  %% ---------- total noise-floor (=harm. mean) in search-band of ringdown frequencies ----------
  SinvSum = zeros ( size ( psd{1}.Sn ) );
  for X = 1:Ndet
    SinvSum += 1 ./ psd{X}.Sn;
    dt{X} = mean( diff ( uvar.ts{X}.ti ) );
    df{X} = mean ( diff ( psd{X}.fk ) );
    fMin{X} = min ( psd{X}.fk );
  endfor
  StotInv = (1/Ndet) * SinvSum;
  Stot = 1./ StotInv;

  assert ( max ( abs ( diff ( [dt{:}] ) )) < 1e-6 );
  assert ( max ( abs ( diff ( [df{:}] ) )) < 1e-6 );
  assert ( max ( abs ( diff ( [fMin{:}] ) )) < 1e-6 );
  dt = dt{1};
  df = df{1};
  fMin = fMin{1};

  Tmax = 5 * max(uvar.prior_tauRange);	%% max time range considered = 5 * tauMax

  %% ---------- templated search over {f0, tau} space ----------
  f0 = [min(uvar.prior_f0Range): uvar.step_f0 : max(uvar.prior_f0Range)];
  tau = [min(uvar.prior_tauRange): uvar.step_tau : max(uvar.prior_tauRange)]';
  [ff0, ttau] = meshgrid ( f0, tau );
  lap_s = 1./ ttau + I * 2*pi * ff0;	%% laplace 'frequency'
  Ntempl = length ( lap_s(:) );
  DebugPrintf (1, "Total number of (f0, tau) templates: %d\n", Ntempl );

  %% ---------- compute 'M_xy = <h_x|h_y>' matrix for SNR term <s|s> ----------
  DebugPrintf ( 1, "Computing M_xy matrix using frequency-domain integral ... ");
  Mxy = compute_Mxy ( fk, ttau, ff0, Stot, Ndet );
  DebugPrintf (1, "done.\n");

  %% ---------- prepare time-shifted and antenna-pattern 'flipped' time-series
  ts = uvar.ts;
  assert ( strcmp ( ts{1}.IFO, "H1" ), "First detector must be 'H1', got '%s'\n", ts{1}.IFO );	%% FIXME, nasty
  assert ( strcmp ( ts{2}.IFO, "L1" ), "Second detector must be 'L1', got '%s'\n", ts{2}.IFO );	%% FIXME, nasty

  shiftBinsL1 = round ( uvar.shiftL1 / dt );
  shiftL1_eff = shiftBinsL1 * dt;
  assert ( abs(shiftL1_eff - uvar.shiftL1) < 1e-6 );
  ts{2}.xi   = - circshift ( ts{2}.xi,   [0, shiftBinsL1] );
  ts{2}.xiW  = - circshift ( ts{2}.xiW,  [0, shiftBinsL1] );
  ts{2}.xiOW = - circshift ( ts{2}.xiOW, [0, shiftBinsL1] );

  yiOW = ts{1}.xiOW + ts{2}.xiOW;

  %% ---------- loop over N start-times [input vector 't0V'] ----------
  numSearches = length ( uvar.t0V );
  clear ("resV");
  for l = 1 : numSearches
    t0_l = uvar.t0V(l);

    %% ----- prepare time-series stresch to analyze
    inds_MaxRange = find ( (ts{1}.ti >= (t0_l - ts{1}.epoch) ) & (ts{1}.ti <= (t0_l - ts{1}.epoch + Tmax)) );
    Dt_i = ts{1}.ti ( inds_MaxRange ) - (t0_l - ts{1}.epoch);
    assert ( min(Dt_i) >= 0 );

    yiOW_l = yiOW ( inds_MaxRange );	%% truncated summed-OW time-series for faster matching

    %% ---------- search parameter-space in {f0, tau} and compute matched-filter in each template ----------
    DebugPrintf ( 1, "Computing match with ringdown templates ... " );
    match = zeros ( size ( lap_s ) );
    for k = 1 : Ntempl	%% loop over all templates
      %% ----- (complex) time-domain template
      hExp_i = exp ( - Dt_i * lap_s(k) );
      match(k) = 2 * dt * sum ( yiOW_l .* hExp_i );
    endfor %% k = 1 : Ntempl
    DebugPrintf ( 1, "done.\n");

    DebugPrintf ( 1, "Computing BSG ... " );
    BSG_f0_tau = zeros ( size ( match ) );
    post_H = zeros ( 1, numH );
    %% marginalize over unknown H scale
    for i = 1 : numH
      H_i      = priorH(i,1);
      priorH_i = priorH(i,2);
      DebugPrintf ( 1, "H = %.2g ...", H_i );
      [ BSG_f0_tau_H, SNR_H{i}, A_H{i}, phi0_H{i} ] = compute_BSG_SNR ( H_i, match, Mxy );
      BSG_f0_tau += priorH_i * BSG_f0_tau_H;	%% marginalize BSG(x;f0,tau,H) over H to get BSG(x;f0,tau)
      post_H(i) = mean ( BSG_f0_tau_H(:) ); 	%% marginalize BSG(x;f0,tau,H) over {f0,tau} to get propto P(H|x)
    endfor
    post_H /= sum (post_H(:));			%% normalize posterior to be sure
    BSG = mean ( BSG_f0_tau(:) );			%% marginalize BSG(x;f0,tau) over {f0,tau} with uniform prior --> BSG(x)
    DebugPrintf ( 1, "done.\n");

    %% pick amplitude-estimates and SNR from H_MPE
    [ val, iMPE ] = max ( post_H(:) );
    H_MPE    = priorH(iMPE,1);
    SNR_est  = SNR_H{iMPE};
    A_est    = A_H{iMPE};
    phi0_est = phi0_H{iMPE};

    bname_l = sprintf ( "Ringdown-GPS%.6fs-f%.0fHz-%.0fHz-tau%.1fms-%.1fms-H%s",
                        t0_l, min(uvar.prior_f0Range), max(uvar.prior_f0Range),
                        1e3 * min(uvar.prior_tauRange), 1e3 * max(uvar.prior_tauRange),
                        priorHname
                      );

    resV(l) = struct ( "bname", bname_l, ...
                       "t0", t0_l, ...
                       "A_est", A_est, ...
                       "phi0_est", phi0_est, ...
                       "BSG", BSG, ...
                       "BSG_f0_tau", {BSG_f0_tau}, ...
                       "SNR", {SNR_est}, ...
                       "post_H", {post_H}, ...
                       "H_MPE", H_MPE ...
                     );

  endfor %% l = 1 : numSearches

  resCommon.Mxy = Mxy;
  resCommon.bname = sprintf ( "Ringdown-f%.0fHz-%.0fHz-tau%.1fms-%.1fms-H%s",
                              min(uvar.prior_f0Range), max(uvar.prior_f0Range),
                              1e3 * min(uvar.prior_tauRange), 1e3 * max(uvar.prior_tauRange),
                              priorHname
                            );
  resCommon.ff0  = ff0;
  resCommon.ttau = ttau;
  resCommon.ts   = ts;
  return

endfunction %% searchRingdown()

function [ BSG, SNR_est, A_est, phi0_est ] = compute_BSG_SNR ( H, match, Mxy )

  Hm2 = H^(-2);
  x_s = - imag ( match );
  x_c =   real ( match );

  det_gamInv = ( Mxy.ss + Hm2 ) .* ( Mxy.cc + Hm2 ) - Mxy.sc.^2;
  det_gam = 1./ det_gamInv;
  normBSG = sqrt(det_gam) * Hm2;

  gam_ss  = det_gam .* ( Mxy.cc + Hm2 );
  gam_cc  = det_gam .* ( Mxy.ss + Hm2 );
  gam_sc  = det_gam .* ( -Mxy.sc );
  x_gam_x = gam_ss .* x_s.^2 + 2 * gam_sc .* x_s .* x_c + gam_cc .* x_c.^2;

  BSG = normBSG .* exp ( 0.5 .* x_gam_x );

  %% ---------- estimate SNR for all templates ----------
  A_s_est = gam_ss .* x_s + gam_sc .* x_c;
  A_c_est = gam_sc .* x_s + gam_cc .* x_c;

  A_est = sqrt ( A_s_est.^2 + A_c_est.^2 );
  phi0_est = - atan2 ( A_s_est, A_c_est );

  SNR2_est = A_s_est.^2 .* Mxy.ss + 2 * A_s_est .* Mxy.sc .* A_c_est + A_c_est.^2 .* Mxy.cc;
  SNR_est = sqrt ( SNR2_est );

  return;
endfunction %% compute_BSG_SNR()
