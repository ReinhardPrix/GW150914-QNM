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

function Mxy = compute_Mxy ( ff0, ttau, psd )

  Ndet = length(psd);
  assert ( size(ttau) == size ( ff0 ) );

  %% ---------- handle buffering of Mxy ----------
  persistent buffer;
  canReuse = false;
  if ( !isempty ( buffer ) )
    same_size_ttau = all((size(ttau) == size(buffer.ttau))(:));
    same_size_ff0  = all((size(ff0) == size(buffer.ff0))(:));
    same_Ndet = (Ndet == length(buffer.psd));
    if ( same_size_ttau && same_size_ff0 && same_Ndet )
      same_ttau = all((ttau == buffer.ttau)(:));
      same_ff0  = all((ff0 == buffer.ff0)(:));
      if (  same_ttau && same_ff0 );
        trueFor = 0;
        for X = 1 : Ndet
          same_IFO = strcmp ( psd{X}.IFO, buffer.psd{X}.IFO );
          same_epoch = (psd{X}.epoch == buffer.psd{X}.epoch);
          same_size_PSD = all ( (size(psd{X}.fk) == size(buffer.psd{X}.fk))(:));
          if ( same_IFO && same_epoch && same_size_PSD )
            same_fk = all ( (psd{X}.fk == buffer.psd{X}.fk)(:) );
            same_Sn = all ( (psd{X}.Sn == buffer.psd{X}.Sn)(:) );
            if (  same_fk && same_Sn )
              trueFor ++;
            endif %% if equal psd{X}
          endif %% if same IFO & epoch & size-PSD
        endfor %% for X = 1:Ndet
        if ( trueFor == Ndet )
          canReuse = true;
        endif %% if true for all X
      endif %% same ttau & ff0
    endif %% same size-ttau&ff0
  endif %% if have buffer

  if ( canReuse )
    DebugPrintf ( 1, "[CAN re-use previous Mxy] ");
    Mxy = buffer.Mxy;
    return;
  else
    DebugPrintf ( 1, "[can NOT re-use a previous Mxy] " );
  endif
  %% ----------

  %% ---------- total noise-floor (=harm. mean) ----------
  fk = psd{1}.fk;
  df = mean ( diff ( fk ) );
  SinvSum = zeros ( size ( psd{1}.Sn ) );
  for X = 1:Ndet
    SinvSum += 1 ./ psd{X}.Sn;
  endfor
  StotInv = (1/Ndet) * SinvSum;
  Stot = 1./ StotInv;

  %% ---------- compute M_xy
  Ntempl = length ( ttau(:) );
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

  %% store results in buffer for potential re-use
  buffer = struct();
  buffer.ff0 = ff0;
  buffer.ttau = ttau;
  buffer.psd = psd;
  buffer.Mxy = Mxy;

  return;

endfunction %% compute_Mxy()
