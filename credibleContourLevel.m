function zConf = credibleContourLevel ( posterior2D, confidence = 0.9 )
  %% return 2D contour iso-pdf level containing 'confidence' probability of the total

  assert  ( ismatrix ( posterior2D ) );
  assert ( (confidence > 0) && (confidence <1 ) );

  norm = sum ( posterior2D (:) );
  posterior2D ./= norm;
  [zConf0, delta, INFO, OUTPUT] = fzero ( @(zIso)  sum ( posterior2D ( find ( posterior2D >= zIso ) )(:) ) - confidence, ...
                                          [ min(posterior2D(:)), max(posterior2D(:)) ], ...
                                          optimset ( "TolX", 1e-4 )
                                        );
  try
    assert ( INFO == 1 );
  catch
    delta
    INFO
    OUTPUT
    error ("fzero() failed\n");
  end_try_catch

  zConf = norm * zConf0;
  return;

endfunction
