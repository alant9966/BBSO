; /**
;  * linear_range.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and plots nonlinearity distribution histograms for all pixels at several potential linear ranges.
;  * Histograms will be saved as .jpg image files in a chosen directory. The console will read "Finished." when the process is finished.
;  * The user may use the histograms to determine what the ideal linear range is to create the nonlinearity curve.
;  * Then, run the 'nonlinearity.pro' file on the data to plot the nonlinearity curve.
;  *
;  * This code assumes that:
;  *    1. The darks and images are taken separately, and each in order of increasing exposure 
;  *        i.e. an example directory would be (1.5 ms img, 3 ms img, 4.5 ms img, 1.5 ms dark, 3 ms dark, 4.5 ms dark) from top to bottom
;  *    2. A uniform light source is used, such that every pixel on the camera is equally and constantly illuminated.
;  *
;  * If one or both of these criteria are not met (if the sun is used instead of uniform illumination), use the "niris_nonlin.pro" code instead,
;  * which ignores irrelevant pixels (outside of the focus area) and accounts for scrambled data (images not in increasing exposure order).      
; */

; Select the directory with the image FITS files
dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

; Derive aperture dimensions from the first FITS file
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

; Create an array mapping every single pixel to its signal over time
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
  img0 = img[*,*, FLOOR((size(img))[3]/2)]
  dk0 = dk[*,*, FLOOR((size(dk))[3]/2)]
  res = fltarr(xsize, ysize)
  
  ; dk0 is not subtracted directly from img0 (instead, each pixel is subtracted individually) to account for
  ;   instances where the pixel's signal may be higher in the dark than the image.
  ; Due to the way IDL works, if the resulting subtracted signal is negative, the value is 'looped' to where saturation
  ;   is considered '0' (e.g. -5 ADU is converted to 65531 ADU). Thus, this code corrects for this issue and sets any
  ;   negative signal values to 0.
  for y = 0, ysize-1 do begin
    for x = 0, xsize-1 do begin
      if (FLOAT(img0[x,y]) - FLOAT(dk0[x,y])) lt 0 then begin
        res[x,y] = 0
      endif else begin
        res[x,y] = FLOAT(img0[x,y] - dk0[x,y])
      endelse
    endfor
  endfor

  ; Add each pixel's signal for the current image to the signals array
  curr = 0UL;
  for y = 0, ysize-1 do begin
    for x = 0, xsize-1 do begin
      signals[i, curr] = res[x,y]
      curr++
    endfor
  endfor
  
  ; The recorded exposure for each image is derived from the header file and added to the exposures array.
  ; Depending on the camera, the index for header[] may need to be changed, as the exposure may be listed in a slightly
  ;   different position in the header file. To do this, view the header file - each parameter listed is counted as one
  ;   index in the array. Count from the top until the index of the exposure parameter is determined (start from index 0).
  exp = header[26] & expsub = strmid(exp, 0, strpos(exp, '.'))
  start = strpos(expsub, ' ', /REVERSE_SEARCH) + 1 
  count = strpos(exp, '/') - start - 1
  exposures[i] = float(strmid(exp, start, count))
endfor

;stop

;=============================== PLOT DISTRIBUTION HISTOGRAMS ===============================

print, 'Plotting nonlinearity histograms...'
print, 'Please pick a directory to save distributions: '
newdir = dialog_pickfile(/directory)

; Plots nonlinearity distributions for half to full range exposures in given intervals.
; For example, if saturation is reached at 2 ms, and the desired interval is 0.1 ms, nonlinearity distributions are plotted for
;   1.0 ms, 1.1 ms, 1.2 ms, 1.3 ms, ... , 2.0 ms. The user can then determine which one is the ideal linear range.
READ, inc, PROMPT='Please provide the desired exposure interval between each nonlinearity histogram (ms). This should be a multiple of the exposure interval: '
interval = exposures[1] - exposures[0]

; For each potential linear range...
for i = FLOOR(n/4)-1, (n/2)-1, FLOOR(inc/interval) do begin
  ; For each pixel...
  for j = 0, (xsize*ysize)-1 do begin
    ; Calculate and add the nonlinearity of the pixel to the nonlins array
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
    ; attempt to divide by 0 and result in -NaN. Instead, immediately mark with a 100% nonlinearity.
    if maxsig eq 0 then begin
      nonlin = 100
    endif else begin
      nonlin = ((maxi+abs(mini))/maxsig)*100
    endelse   
    nonlins[j] = nonlin
  endfor  
  
  ; To keep the histograms clean and readable to the user, treat extreme outliers (nonlinearity > 50%) as missing data.
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