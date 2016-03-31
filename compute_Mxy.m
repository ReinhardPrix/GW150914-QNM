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

function Mxy = compute_Mxy ( fk, ttau, ff0, Stot, Ndet )

  assert ( length ( fk ) == length ( Stot ) );
  assert ( size(ttau) == size ( ff0 ) );
  Ntempl = length ( ttau(:) );
  df = mean ( diff ( fk ) );

  M_ss = M_cc = M_sc = zeros ( size ( ttau ) );
  for l = 1 : Ntempl
    %% ----- whitened frequency-domain template basis functions ----------
    denom_k = ( 1 + I * 4*pi * fk * ttau(l) - 4*pi^2 * ( fk.^2 - ff0(l)^2 ) * ttau(l)^2 );
    hsFT_k  = ttau(l) * ( 2*pi * ff0(l) * ttau(l) ) ./ denom_k;
    hcFT_k  = ttau(l) * ( 1 + I * 2*pi*fk * ttau(l) ) ./ denom_k;
    %% ----- compute M-matrix from template self-match integrals in frequency-domain ----------
    M_ss(l) = 4 * Ndet * df * sum ( abs(hsFT_k).^2 ./ Stot );
    M_cc(l) = 4 * Ndet * df * sum ( abs(hcFT_k).^2 ./ Stot );
    M_sc(l) = 4 * Ndet * df * real ( sum ( hsFT_k .* conj(hcFT_k) ./ Stot ) );
  endfor %% l

  Mxy.ss = M_ss;
  Mxy.cc = M_cc;
  Mxy.sc = M_sc;
  return;

endfunction %% compute_Mxy()
