; /**
;  * mask.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and identifies and outputs the coordinates of hot pixels and dead pixels.
;  * Coordinates of defect pixels are also stored in corresponding arrays for referencing.
;  *
;  * This code assumes that the darks and images are taken separately, and each in order of increasing exposure
;  *    i.e. an example directory would be (1.5 ms img, 3 ms img, 4.5 ms img, 1.5 ms dark, 3 ms dark, 4.5 ms dark)
;  *
;  * If this criteria is not met, use the "niris_mask" code instead, which accounts for scrambled data (images not in increasing exposure order).
; */

; Select the directory with the image FITS files
dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

; Derive image dimensions from the first FITS file
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
signals = fltarr(n/2, (xsize*ysize)) & imgsig = signals
mask = fltarr(xsize, ysize)

; Arrays holding the coordinates of nonlinear, dead, and hot pixels
nonlins = fltarr(xsize, ysize) & dead = nonlins & hot = nonlins

;=============================== READ AND ANALYZE IMAGES ===============================

; For each image...
for i = 0, (n/2)-1 do begin
  ; Subtract dark images, choosing the middle frame
  img = readfits(list[i], header)
  dk = readfits(list[i+(n/2)])

  img0 = img[*,*, FLOOR((size(img))[3]/2)]
  dk0 = dk[*,*, FLOOR((size(dk))[3]/2)]
  
  ; dk0 is not subtracted directly from img0 (instead, each pixel is subtracted individually) to account for
  ;   instances where the pixel's signal may be higher in the dark than the image.
  ; Due to the way IDL works, if the resulting subtracted signal is negative, the value is 'looped' to where saturation
  ;   is considered '0' (e.g. -5 ADU is converted to 65531 ADU). Thus, this code corrects for this issue and sets any
  ;   negative signal values to 0.
  res = fltarr(xsize, ysize)
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
      ; Post-dark subtraction
      signals[i, curr] = res[x,y]
      ; Pre-dark subtraction
      imgsig[i, curr] = img0[x,y]
      curr++
    endfor
  endfor

  ; The recorded exposure for each image is added to the exposures array
  exp = header[26] & expsub = strmid(exp, 0, strpos(exp, '.'))
  start = strpos(expsub, ' ', /REVERSE_SEARCH) + 1
  count = strpos(exp, '/') - start - 1
  exposures[i] = float(strmid(exp, start, count))
endfor

READ, linrange, PROMPT='Please enter the linear range (ADU): '

linindex = (n/2)-1
; Determines the index to regress over (the linear range) based on the inputted value.
; The average signal atn each exposure time is manually calculated to ignore pixels outside of the focus area (if the sun is used for testing).
; If a uniform light source is used, this functions the exact same as the avg() method.
for i = 0, (n/2)-2 do begin
  avgsig = 0.0
  avgsignext = 0.0
  count = 0.0
  
  ; Only uses pixels that saturate or nearly saturate over the entire testing range (all pixels that are illuminated) 
  ;   to be used in calculating the average.
  ; Here, I've determined that any pixels that increase more than 10000 ADU between the lowest and highest exposures tested
  ;   must be pixels within the area of interest (pixels outside of the viewing area of the sun stay at near-zero signal at 
  ;   any exposure), as the saturation value is 65536 ADU.
  ; This threshold may need to be changed for cameras with a different bit-depth and saturation value. For example, for a 12-bit
  ;   camera may want a threshold of 400 ADU instead. This is not majorly important as all pixels of interest should fully saturate
  ;   while all other pixels should always stay at nearly zero signal, so the threshold only needs to be below the saturation value.
  for k = 0, xsize*ysize-1 do begin
    if ~(avg(signals[*,k]) lt 10000) then begin
      avgsig += signals[i,k]
      avgsignext += signals[i+1,k]
      count++
    endif
  endfor
  avgsig /= count
  avgsignext /= count
  
  ; At each index, checks if the inputted linear range (in ADU) is between the average signal of the current index and
  ;   the average signal of the next index.
  if (avgsig le linrange) and (avgsignext ge linrange) then begin
    linindex = i
    break
  endif
  print, i
endfor

interval = exposures[1] - exposures[0]

;====================== NONLINEARITY ======================

print, 'Coordinates of nonlinear pixels: '

; For each pixel, calculate the nonlinearity and determine if the pixel is exceedingly nonlinear.
for i = 0, (xsize*ysize)-1 do begin
  ; Perform a linear regression and calculate nonlinearity
  coeff = LINFIT(exposures[0:linindex], signals[0:linindex, i])
  mini = signals[linindex,i]
  maxi = 0.0
  maxsig = 0.0
  for index = 0, linindex do begin
    curr = signals[index,i] - (coeff[1]*exposures[index] + coeff[0])
    if curr gt maxi then begin
      maxi = curr
    endif
    if curr lt mini then begin
      mini = curr
    endif
    if signals[index,i] gt maxsig then begin
      maxsig = signals[index,i]
    endif
  endfor

  ; If a pixel is completely dead (max signal of 0), do not calculate nonlinearity as it will
  ; attempt to divide by 0 and result in NaN. Instead, immediately mark with a 100% nonlinearity.
  if maxsig eq 0 then begin
    nonlin = 100
  endif else begin
    nonlin = ((maxi+abs(mini))/maxsig)*100
  endelse

  ; Determine the (x,y) coordinates of the current pixel
  xc = i MOD xsize
  yc = FLOOR(i/xsize)

  ; If the pixel has a nonlinearity greater than 10%, it is recorded as a particularly nonlinear pixel and printed to the console.
  ; Depending on what you consider to be 'overly nonlinear', this threshold value is not concrete, and should be changed.
  if (nonlin gt 10) then begin
    mask[xc, yc] = 1
    nonlins[xc, yc] = nonlin
    print, '('+strtrim(string(xc),2)+', '+strtrim(string(yc),2)+')
  endif else begin
    mask[xc, yc] = 0
  endelse
endfor

;====================== DEAD PIXELS AND HOT PIXELS ======================

print, 'Coordinates of dead pixels and hot pixels: '

; Testing for hot pixels - checks if the pixel saturates unnaturally quickly.
; Testing for dead pixels - checks if the pixel is nearly zero when it should be saturated.
hotindex = FLOOR(linindex/2)
for i = 0, (xsize*ysize)-1 do begin
  ; (x,y) coordinates
  xc = i MOD xsize
  yc = FLOOR(i/xsize)

  ; At the highest exposure tested, the current pixel should be saturated. If not, and if the signal is sufficiently near zero, 
  ;   mark it as a dead pixel and print its coordinates to the console.
  ; Depending on what you consider to be 'near zero', this threshold value is not concrete, and should be changed.
  ;   However, this threshold value must not be greater than the saturation value.
  if imgsig[(n/2)-1,i] lt 10000 then begin
    dead[xc,yc] = 1
    mask[xc,yc] = 1
    print, 'dead: ('+strtrim(string(xc),2)+', '+strtrim(string(yc),2)+')
  endif

  ; At the index of half the linear range, the pixel should not be saturated (it should show a signal of roughly saturation/2).
  ;   If the pixel is saturated at this index, mark it as a hot pixel and print its coordinates to the console.
  ; Depending on your definition of a hot pixel, the hotindex variable should be modified to look for pixels that saturate even quicker, or slower.
  if imgsig[hotindex,i] ge avg(signals[(n/2)-1,i]) then begin
    hot[xc,yc] = 1
    mask[xc,yc] = 1
    print, 'hot: ('+strtrim(string(xc),2)+', '+strtrim(string(yc),2)+')
  endif
endfor

end