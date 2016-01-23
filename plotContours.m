function plotContours ( in )
  global iFig0 = 0;
  global tMergerOffs
  global tEvent;
  global f0GR;
  global taumsGR;

  Nsteps = length ( in );

  figure ( iFig0 + 5 ); clf;
  hold on;
  he = errorbar ( f0GR.val, taumsGR.val, f0GR.lerr, f0GR.uerr, taumsGR.lerr, taumsGR.uerr, "~>.r;IMR;" );
  for i = 1 : Nsteps
    [C, H] = contour ( in{i}.ff0, in{i}.ttau * 1e3, in{i}.posterior2D, in{i}.isoConf2 * [ 1, 1 ] );
    set ( H, "linecolor", "black", "linewidth", 2 );
    leg = sprintf ( "tM + %.1fms", (in{i}.tGPS - tEvent - tMergerOffs) * 1e3 );
    cl = clabel ( C, H, "FontSize", 12, "Color", 'k');
    set ( cl, "string", "" );
    %%set ( cl(1), "string", leg );
    set ( cl(end), "string", leg );
  endfor
  xlabel ( "Freq [Hz]" );
  ylabel ( "tau [ms]" );
  hold off;
  grid on;

  fname = sprintf ( "%s-contours.pdf", in{1}.bname );
  ezprint ( fname, "width", 512 );

  return;

endfunction
