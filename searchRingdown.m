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

  uvar = parseOptions ( varargin,
                        {"multiTS", "cell" },	%% cell-array [over detectors]: normal, whitentend, and over-whitened timeseries
                        {"multiPSD", "cell"},	%% cell-array [over detectors]: PSD estimate over frequency range, for each detector
                        {"t0V", "real,strictpos,vector", 1126259462.43 },	%% ringdown start-time in GPS s
                        {"prior_f0Range", "real,strictpos,vector", [220,  270] },
                        {"step_f0", "real,strictpos,scalar", 0.1 },
                        {"prior_tauRange", "real,vector", [1e-3, 30e-3] },
                        {"step_tau", "real,strictpos,scalar", 0.2e-3 },
                        {"prior_H", "real,strictpos,matrix", [4e-22, 1]},
                        {"plotResults", "bool", false },
                        {"skyCorr", "struct"}		%% detector-dependent time-shift and antenna-pattern infos
                      );

  numIFOs = length ( uvar.multiTS );
  assert ( (numIFOs == length(uvar.multiPSD)) && (numIFOs == length(uvar.skyCorr)) );

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

  dt = 1 / uvar.multiTS{1}.fSamp;
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
  Mxy = compute_Mxy ( ff0, ttau, uvar.multiPSD );
  DebugPrintf (1, "done.\n");

  %% ---------- prepare time-shifted and antenna-pattern corrected *summed* time-series
  tsSum.yiOW = zeros ( size ( uvar.multiTS{1}.xi ) );
  multiTScorr = uvar.multiTS;
  for X = 1 : numIFOs
    IFO = uvar.multiTS{X}.IFO
    X1 = find ( strcmp ( IFO, {uvar.skyCorr.IFO} ) ); assert ( !isempty(X1) && (length(X1) == 1) );
    shiftBins = round ( uvar.skyCorr(X1).timeShift / dt );
    ampFact   = uvar.skyCorr(X1).ampFact;
    shift_eff = shiftBins * dt; assert ( abs(shift_eff - uvar.skyCorr(X1).timeShift) < 1e-6 );	%% check we're within 0.001ms

    %% there's 2 ways to do the time-shifts: add shift-value to ti, or shift bins in xi
    %% for now we use the second method, as we're going to sum the time-shifted ts{X} together and
    %% compute a single Bayes-factor.
    multiTScorr{X}.xi   = ampFact * circshift ( uvar.multiTS{X}.xi,   [0, shiftBins] );
    multiTScorr{X}.xiW  = ampFact * circshift ( uvar.multiTS{X}.xiW,  [0, shiftBins] );
    multiTScorr{X}.xiOW = ampFact * circshift ( uvar.multiTS{X}.xiOW, [0, shiftBins] );

    tsSum.yiOW += multiTScorr{X}.xiOW;
  endfor %% X = 1 : Net
  tsSum.ti = multiTScorr{1}.ti;
  tsSum.epoch = multiTScorr{1}.epoch;

  %% ---------- loop over N start-times [input vector 't0V'] ----------
  numSearches = length ( uvar.t0V );
  clear ("resV");
  for l = 1 : numSearches
    t0_l = uvar.t0V(l);

    %% ----- prepare time-series stresch to analyze
    inds_MaxRange = find ( (tsSum.ti >= (t0_l - tsSum.epoch) ) & (tsSum.ti <= (t0_l - tsSum.epoch + Tmax)) );
    Dt_i = tsSum.ti ( inds_MaxRange ) - (t0_l - tsSum.epoch);
    assert ( min(Dt_i) >= 0 );

    yiOW_l = tsSum.yiOW ( inds_MaxRange );	%% truncated summed-OW time-series for faster matching
    %% ---------- search parameter-space in {f0, tau} and compute matched-filter in each template ----------
    DebugPrintf ( 1, "----- t0 = %.6f s ----------\n", t0_l );
    DebugPrintf ( 1, "Computing match with ringdown templates ... " );
    match = zeros ( size ( lap_s ) );
    for k = 1 : Ntempl	%% loop over all templates
      %% ----- (complex) time-domain template
      hExp_i = exp ( - Dt_i * lap_s(k) );
      match(k) = 2 * dt * sum ( yiOW_l .* hExp_i );
    endfor %% k = 1 : Ntempl
    DebugPrintf ( 1, "done.\n");

    DebugPrintf ( 1, "Computing BSG(f0,tau) ... " );
    BSG_f0_tau = zeros ( size ( match ) );
    post_H = zeros ( 1, numH );
    %% marginalize over unknown H scale
    for i = 1 : numH
      H_i      = priorH(i,1);
      priorH_i = priorH(i,2);
      [ BSG_f0_tau_H, AmpV_ML{i} ] = compute_BSG ( H_i, match, Mxy );
      BSG_f0_tau += priorH_i * BSG_f0_tau_H;	%% marginalize BSG(x;f0,tau,H) over H to get BSG(x;f0,tau)
      post_H(i) = mean ( BSG_f0_tau_H(:) ); 	%% marginalize BSG(x;f0,tau,H) over {f0,tau} to get propto P(H|x)
    endfor %% i = 1:numH
    post_H /= sum (post_H(:));			%% normalize posterior to be sure
    BSG = mean ( BSG_f0_tau(:) );			%% marginalize BSG(x;f0,tau) over {f0,tau} with uniform prior --> BSG(x)
    DebugPrintf ( 1, "done.\n");

    [ BSG_MP, iMP ] = max ( BSG_f0_tau(:) );
    f0MP = ff0(iMP);
    tauMP = ttau(iMP);
    %% for plotting OW-timeseries re-scaled to MP template: store noise-values at MP frequency
    SX_MP = zeros ( 1, numIFOs );
    for X = 1 : numIFOs
      [val, freqInd] = min ( abs ( uvar.multiPSD{X}.fk - f0MP ) );
      SX_MP(X) = uvar.multiPSD{X}.Sn ( freqInd );
      DebugPrintf ( 2, "X = %d: sqrt(SX)|_MP = %g /sqrt(Hz)\n", X, sqrt ( SX_MP(X) ) );
    endfor
    lambdaMP = struct ( "iMP", iMP, "BSG_MP", BSG_MP, "f0", f0MP, "tau", tauMP, "SX", SX_MP );
    match_k = match(iMP);
    Mxy_k = struct ( "ss", Mxy.ss(iMP), "cc", Mxy.cc(iMP), "sc", Mxy.sc(iMP) );

    %% ---------- obsolete: SNR from naive maximum-likelihood amplitude estimate ----------
    [ val, iH_MPE ] = max ( post_H(:) );	%% pick amplitude-estimates and SNR from H_MPE
    A_s_ML  = AmpV_ML{iH_MPE}.A_s(iMP);
    A_c_ML  = AmpV_ML{iH_MPE}.A_c(iMP);
    A_ML    = sqrt ( A_s_ML^2 + A_c_ML^2 );
    phi0_ML = - atan2 ( A_s_ML, A_c_ML );
    SNR_ML  = sqrt ( A_s_ML^2 * Mxy_k.ss + 2 * A_s_ML * Mxy_k.sc * A_c_ML + A_c_ML^2 * Mxy_k.cc );
    AmpML = struct ( "A_s", A_s_ML, "A_c", A_c_ML, "A", A_ML, "phi0", phi0_ML, "SNR", SNR_ML );

    %% ---------- maximum-posterior amplitude estimate in MPE-point {f0,tau} ----------
    DebugPrintf ( 1, "Estimating MP amplitudes (A,phi0) and SNR in MP(f0,tau) ... ");
    maxH = max ( priorH(:,1) );
    vA = linspace ( -maxH, maxH, 100 );
    [A_s, A_c] = meshgrid ( vA, vA );
    post_vA = zeros ( size ( A_s ) );
    for i = 1 : numH
      H_i      = priorH(i,1);
      priorH_i = priorH(i,2);
      post_vA_H = posterior_A_at_lambda_H ( A_s, A_c, H_i, match_k, Mxy_k );
      post_vA += priorH_i * post_vA_H;
    endfor %% i = 1:numH
    [maxPA, iMaxPA] = max ( post_vA(:) );
    A_s_MP  = A_s(iMaxPA);
    A_c_MP  = A_c(iMaxPA);
    A_MP    = sqrt ( A_s_MP^2 + A_c_MP^2 );
    phi0_MP = - atan2 ( A_s_MP, A_c_MP );
    SNR_MP  = sqrt ( A_s_MP^2 * Mxy_k.ss + 2 * A_s_MP * Mxy_k.sc * A_c_MP + A_c_MP^2 * Mxy_k.cc );

    AmpMP = struct ( "A_s", A_s_MP, "A_c", A_c_MP, "A", A_MP, "phi0", phi0_MP, "SNR", SNR_MP );
    DebugPrintf ( 1, "done.\n");

    bname_l = sprintf ( "Ringdown-GPS%.6fs-f%.0fHz-%.0fHz-tau%.1fms-%.1fms-H%s",
                        t0_l, min(uvar.prior_f0Range), max(uvar.prior_f0Range),
                        1e3 * min(uvar.prior_tauRange), 1e3 * max(uvar.prior_tauRange),
                        priorHname
                      );

    resV(l) = struct ( "bname", bname_l, ...
                       "t0", t0_l, ...
                       "BSG", BSG, ...
                       "BSG_f0_tau", {BSG_f0_tau}, ...
                       "post_H", {post_H}, ...
                       "match", {match}, ...
                       "lambdaMP", lambdaMP, ...
                       "AmpMP", AmpMP, ...
                       "AmpML", AmpML ...
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
  resCommon.multiTScorr = multiTScorr;
  return

endfunction %% searchRingdown()

function [ BSG, AmpML ] = compute_BSG ( H, match, Mxy )

  Hm2 = H^(-2);
  x_s = - imag ( match );
  x_c =   real ( match );

  det_gamInv = ( Mxy.ss + Hm2 ) .* ( Mxy.cc + Hm2 ) - Mxy.sc.^2;
  det_gam = 1 ./ det_gamInv;
  normBSG = sqrt ( det_gam ) * Hm2;

  gam_ss  = det_gam .* ( Mxy.cc + Hm2 );
  gam_cc  = det_gam .* ( Mxy.ss + Hm2 );
  gam_sc  = det_gam .* ( -Mxy.sc );
  x_gam_x = gam_ss .* x_s.^2 + 2 * gam_sc .* x_s .* x_c + gam_cc .* x_c.^2;

  ln_BSG = log(normBSG) + ( 0.5 .* x_gam_x );
  BSG = exp ( ln_BSG );

  %% ---------- simple maximum-likelihood amplitude estimate over all templates ----------
  AmpML.A_s = gam_ss .* x_s + gam_sc .* x_c;
  AmpML.A_c = gam_sc .* x_s + gam_cc .* x_c;

  return;
endfunction %% compute_BSG_SNR()


function post_vA = posterior_A_at_lambda_H ( As, Ac, H, match_k, Mxy_k )
  assert ( size ( As )  == size ( Ac ) );

  Hm2 = H^(-2);
  x_s = - imag ( match_k );
  x_c =   real ( match_k );

  gamInv.ss = Mxy_k.ss + Hm2;
  gamInv.cc = Mxy_k.cc + Hm2;
  gamInv.sc = Mxy_k.sc;

  A_gamInv_A = As .* gamInv.ss .* As + 2 * As .* gamInv.sc .* Ac + Ac .* gamInv.cc .* Ac;
  A_x = As * x_s + Ac * x_c;

  ln_post_vA = - 0.5 * A_gamInv_A + A_x;
  post_vA = exp ( ln_post_vA );
  return;

endfunction %% posterior_A_at_lambda_H
