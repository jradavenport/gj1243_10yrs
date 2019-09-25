pro fluxcal
  set_plot,'X'
  plotstuff,/set,/silent

  ; .r /Users/james/idl/pro/remove
  ; .r /Users/james/idl/pro/planck
  
  loadct,39,/silent
  ;===== compute the flux calibration for this target

  
  ;--> need a flux calibration spectrum: PMSU to the Rescue!

  readcol, 'GJ1243_PMSU1.txt', wave_cal, flux_cal, /silent
  remove,where(wave_cal gt 7500), wave_cal, flux_cal

  ;; sort spectrum, if needed...
  ;ss = sort(wave_cal)
  ;wave_cal = wave_cal[ss]
  ;flux_cal = flux_cal[ss]

  
  ; Bochanski M4 active template
  ;std = mrdfits('m4.all.na.k.fits',0,hdr,/silent)
  ;std_flux = std[*,1]
  ;std_wave = sxpar(hdr, 'CRVAL1') + $
  ;           findgen(n_elements(std_flux)) * sxpar(hdr, 'CD1_1')


  
  ; the Davenport+2012 "oir_template"
  readcol,'m4_template.dat', std2_wave, std2_flux,/silent
  remove, where(std2_flux lt -10), std2_flux, std2_wave
  std_flux = 10^std2_flux
  std_wave = std2_wave * 10000
  
  ; now calibrate the template to this spectrum
  x1 = where(wave_cal ge 7200 and wave_cal le 7450)
  x2 = where(std_wave ge 7200 and std_wave le 7450)
  flux_convert = median(flux_cal[x1]) / median(std_flux[x2])

  
  plot, wave_cal, flux_cal,xrange=[4500,10500]
  oplot, std_wave, std_flux*flux_convert, color=250
  ;oplot, std2_wave, (std2_flux)*flux_convert2, color=150
  
  ; the TESS filter curve
  readcol, 'tess-response-function-v1.0.csv', w_f, trans_f,/silent
  wave_f = w_f * 10.
  oplot, wave_f, trans_f * max(flux_cal), color=90

  ; the Kepler filter curve
  readcol, 'kepler_response_hires1.txt', w_fk, trans_fk,/silent
  wave_fk = w_fk * 10
  trans_fk =  trans_fk/max(trans_fk)
  oplot, wave_fk, trans_fk*max(flux_cal), color=200


  
  ; - clip spectrum to filter wavelength range
  xok = where(std_wave ge min(wave_f) and $
              std_wave le max(wave_f))
  std_wave2 = std_wave[xok]
  std_flux2 = std_flux[xok] * flux_convert

  xok = where(std_wave ge min(wave_fk) and $
              std_wave le max(wave_fk))
  std_wave2k = std_wave[xok]
  std_flux2k = std_flux[xok] * flux_convert


  ;- resample filter on to spectrum wavelength grid
  trans_f2 = interpol(trans_f,wave_f,std_wave2)
  trans_f2k = interpol(trans_fk,wave_fk,std_wave2)

  x50 = where(trans_f2 ge max(trans_f2)/2.)
  FWHM = max(std_wave2[x50]) - min(std_wave2[x50]) ; in Angstroms

  x50k = where(trans_f2k ge max(trans_f2k)/2.)
  FWHMk = max(std_wave2[x50k]) - min(std_wave2[x50k]) ; in Angstroms


  plx = 83.4814  ; mas
  dist = 1000. / plx * !pc

  print,'Distance (pc) = ',dist / !pc

  ;- convolve spectrum with filter
 
  L_POINT = alog10(TSUM(std_wave2, std_flux2 * trans_f2) * $
                   (2. * !dpi * dist^2.))

  L_POINTk = alog10(TSUM(std_wave2k, std_flux2k * trans_f2k) * $
                    (2. * !dpi * dist^2.))

  print,'log LUMINOSITY = ',l_point, l_pointk


  bb10_tess = planck(std_wave2, 10000)
  bb10_kepl = planck(std_wave2k, 10000)

  print, total(bb10_tess) / total(bb10_kepl)


  set_plot,'ps'
  device, filename='gj1243_spectrum.eps',/color,/encap,/inch,xsize=8,ysize=5
  plot, std_wave/10, std_flux*flux_convert/1d-13, $
        xrange=[350,1150],xsty=1, /nodata, $
        xtitle='Wavelength (nm)', ytitle='Flux (10!u-13!n erg s!u-1!ncm!u-2!nA!u-1!n)'
  loadct,0,/silent
  oplot, std_wave/10, std_flux*flux_convert/1d-13, thick=1,color=120
  loadct,39,/silent
  oplot, wave_cal/10, flux_cal/1d-13, thick=8
  oplot, wave_f/10, trans_f , color=75, thick=5
  oplot, wave_fk/10, trans_fk , color=250, thick=5
  device,/close

  set_plot,'X'

  
  return
end
