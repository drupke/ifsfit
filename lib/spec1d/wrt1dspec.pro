;+
; NAME:
;	WRT1DSPEC
;
; PURPOSE:
;	Write an RKSPEC-format one-dimensional spectrum.
;
; INPUTS:
;	specname - output file name
;	spec     - spectrum
;	sigspec  - sigma spectrum
;	sky      - sky spectrum
;	mask     - bad pixel mask
;	header   - FITS header
;
; OPTIONAL INPUTS:
;	datapath - I/O data path
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMENTS:
;
; PROCEDURES USED:
;	MWRFITS, CWD()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 January 21, U of A
;-

pro wrt1dspec, specname, spec, sigspec, sky, mask, $
  header, datapath=datapath

    if not keyword_set(datapath) then datapath = cwd()

    mwrfits, spec, datapath+specname, header, /create
    mwrfits, sigspec, datapath+specname
    mwrfits, sky, datapath+specname
    mwrfits, mask, datapath+specname

return
end
