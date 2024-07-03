; /**
;  * linearity_ALL.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and plots the nonlinearity distribution for all pixels.
;  * A separate IDL program is run to create a mask for nonlinear and defect pixels.
; */

dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

; Derive image dimensions from the first FITS file
print, 'Initializing parameters and variables...'
img0 = readfits(list[0], head)
xlen = head[3] & xsub = strmid(xlen, 0, strpos(xlen, '/')-2)
ylen = head[4] & ysub = strmid(ylen, 0, strpos(ylen, '/')-2)
xstart = strpos(xsub, ' ', /REVERSE_SEARCH) + 1
ystart = strpos(ysub, ' ', /REVERSE_SEARCH) + 1
xcount = strpos(xlen, '/') - xstart - 1
ycount = strpos(ylen, '/') - ystart - 1
xsize = float(strmid(xlen, xstart, xcount))
ysize = float(strmid(ylen, ystart, ycount))

; Create an array mapping every single pixel (~4 million) to its signal over time
; (px1)[2000, 2031, 2403, ...]
; (px2)[2000, 2031, 2403, ...]
; (px3)[2000, 2031, 2403, ...]
; ...
exposures = fltarr(n/2)
signals = fltarr(n/2, (xsize*ysize))
nonlins = fltarr(xsize*ysize)

;=============================== READ AND ANALYZE IMAGES ===============================

print, 'Reading images...'
; For each image...
for i = 0, (n/2)-1 do begin
  ; Subtract dark images, choosing the middle frame
  img = readfits(list[i], header)
  dk = readfits(list[i+(n/2)])
  res = img - dk
  
  ; Add each pixel's signal for the current image to the signals array
  curr = 0UL;
  res0 = res[*,*, FLOOR((size(res))[3]/2)]
  for y = 0, ysize-1 do begin
    for x = 0, xsize-1 do begin
      signals[i, curr] = res0[x,y]
      ;print, 'Inserting ', strtrim(string(res0[x,y])), ' at (', strtrim(string(i)), ',',  strtrim(string(curr)), ')'
      ;wait, 1
      curr++
    endfor
  endfor  
  
  ; The recorded exposure for each image is added to the exposures array.
  exp = header[26] & expsub = strmid(exp, 0, strpos(exp, '.'))
  start = strpos(expsub, ' ', /REVERSE_SEARCH) + 1 
  count = strpos(exp, '/') - start - 1
  exposures[i] = float(strmid(exp, start, count))
endfor

;stop

print, 'Plotting nonlinearity histograms...'
print, 'Please pick a directory to save distributions: '
newdir = dialog_pickfile(/directory)
; Plots nonlinearity distributions for half to full range exposures in increments of 0.1 ms
interval = exposures[1] - exposures[0]
; For each possible linear range...
for i = FLOOR(n/4)-1, (n/2)-1, FLOOR(0.1/interval) do begin
  ; For each pixel...
  for j = 0, (xsize*ysize)-1 do begin
    ; Calculate and add the nonlinearity to the corresponding array
    coeff = LINFIT(exposures[0:i], signals[0:i, j]) 
    mini = signals[i,j]
    maxi = 0.0
    maxsig = 0.0
    for index = 0, i do begin
      curr = signals[index,j] - (coeff[1]*exposures[index] + coeff[0])
      if curr gt maxi then begin
        maxi = curr
      endif
      if curr lt mini then begin
        mini = curr
      endif
      if signals[index,j] gt maxsig then begin
        maxsig = signals[index,j]
      endif
    endfor
    
    ; If a pixel is completely dead (max signal of 0), do not calculate nonlinearity as it will
    ; attempt to divide by 0 and result in NaN. Instead, immediately mark with a 100% nonlinearity.
    if maxsig eq 0 then begin
      nonlin = 100
    endif else begin
      nonlin = ((maxi+abs(mini))/maxsig)*100
    endelse   
    nonlins[j] = nonlin
  endfor  
  
  for ind = 0, (size(nonlins))[1]-1 do begin
    if nonlins[ind] gt 50 then begin
      nonlins[ind] = !VALUES.F_NAN
    endif
  endfor
  
  ; Plot and save the distribution histogram to a specific directory
  avgsig = avg(signals[i,*])
  hist = HISTOGRAM(nonlins, LOCATIONS=binvals, BINSIZE=0.1, /NAN)
  histplot = PLOT(binvals, hist, XTITLE='Nonlinearity (%)', YTITLE='Frequency', $
    TITLE='Pixel Nonlinearity Distributions ('+strtrim(string(avgsig),2)+' ADU)', /STAIRSTEP)
  savedir = newdir+'/nonlin_'+strtrim(string(avgsig),2)+'adu.jpg'
  histplot.Save, savedir
endfor

print, 'Finished.'

end