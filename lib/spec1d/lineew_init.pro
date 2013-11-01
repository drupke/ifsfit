function lineew_init, linewave, linename
; function for RKSPLOTEW()
    
    lineew = create_struct('linename', linename, 'linewave', $
      linewave, 'continuum', [0.0,-1.0], 'boxflux', [0.0,-1.0], $
      'gaussflux', [0.0,-1.0], 'ew', [0.0,-1.0])

return, lineew
end
