function inds = binRange ( xMin, xMax, xVals )
  %% unified handling of finding 'bin ranges' from bounds 'xMin|xMax' on sets of continuous quantities
  %% 'xVals', in such a way to be robust wrt roundoff-errors.
  %% Here we simply round to the closest bin, which is OK when it's not important to guarantee
  %% that we always include the bounds 'xMin|xMax'

  [val, iMin] = min ( abs ( xVals -  xMin ) );
  [val, iMax] = min ( abs ( xVals -  xMax ) );

  inds = [iMin : iMax ];

  return;
endfunction
