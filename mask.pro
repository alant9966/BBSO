; /**
;  * mask.pro
;  * Alan Tong
;  *
;  * Identifies and outputs the coordinates of nonlinear and defect pixels to create a mask.
; */

; PROBLEM: even if pixle is dead or hot, if there is a very minor but consistent linear increase it will pass as a non-defective pixel
; SOLUTION: add in separate standards for dead and hot pixels? maybe if the difference between the starting and ending is too small mark it as defective
    ; ... or perform another slope hypothesis test lol
    
;!EXCEPT = 0

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
signals = fltarr(n/2, (xsize*ysize))
mask = fltarr(xsize, ysize)

; Arrays for troubleshooting/testing purposes:
nonlins = fltarr(xsize, ysize) & pvals = nonlins & incs = nonlins

;=============================== READ AND ANALYZE IMAGES ===============================

; For each image...
for i = 0, (n/2)-1 do begin
  ; Subtract dark images
  img = readfits(list[i], header)
  dk = readfits(list[i+(n/2)])
  res = img - dk

  ; Add each pixel's signal for the current image to the signals array
  curr = 0UL;
  res0 = res[*,*, FLOOR((size(res))[3]/2)]
  for y = 0, ysize-1 do begin
    for x = 0, xsize-1 do begin
      signals[i, curr] = res0[x,y]
      curr++
    endfor
  endfor 

  ; The recorded exposure for each image is added to the exposures array.
  exp = header[26] & expsub = strmid(exp, 0, strpos(exp, '.'))
  start = strpos(expsub, ' ', /REVERSE_SEARCH) + 1
  count = strpos(exp, '/') - start - 1
  exposures[i] = float(strmid(exp, start, count))
endfor

; NOTE: the linear range should be a multiple of the exposure interval
READ, linrange, PROMPT='Please enter the linear range (ADU): '
linindex = 0
for i = 0, (n/2)-2 do begin
  if (avg(signals[i,*]) le linrange) and (avg(signals[i+1,*]) ge linrange) then begin
    linindex = i
    break
  endif
endfor

interval = exposures[1] - exposures[0]

;; For each pixel...
;for i = 0, (xsize*ysize) do begin
;  ; Calculate and add nonlinearities and correlation coefficients to the corresponding arrays
;  coeff = LINFIT(exposures[0:linindex], signals[0:linindex, i])
;  corr = CORRELATE(exposures[0:linindex], signals[0:linindex, i])
;  
;  mini = signals[linindex,i]
;  maxi = 0.0
;  maxsig = 0.0
;  for index = 0, linindex do begin
;    curr = res[index,i] - (coeff[1]*exposures[index] + coeff[0])
;    if curr gt maxi then begin
;      maxi = curr
;    endif
;    if curr lt mini then begin
;      mini = curr
;    endif
;    if signals[index,i] gt maxsig then begin
;      maxsig = signals[index,i]
;    endif
;  endfor
;  
;  nonlin = ((maxi+abs(mini))/maxsig)*100
;  nonlins[i] = nonlin
;  corrs[i] = corr
;endfor

;====================== NONLINEARITY AND CORRELATION-COEFFICIENTS ======================

print, 'Coordinates of defect/nonlinear pixels: '
mx = (RUNNING_STATS(signals[linindex,*]))[0]

; For each pixel, calculate the nonlinearity and correlation coefficient.
; Then, determine if the pixel should be masked and add its coordinates to an array.
for i = 0, (xsize*ysize)-1 do begin
  coeff = LINFIT(exposures[0:linindex], signals[0:linindex, i])
  corr = CORRELATE(exposures[0:linindex], signals[0:linindex, i])

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
  
  ;----- Calculating the p-value of the OLS linear regression using a statistical t-test -----

  n = linindex + 1
  X = fltarr(n,2)
  X[*, 0] = 1.0
  X[*, 1] = exposures[0:linindex]
  
  XT = TRANSPOSE(X)
  XTX = XT # X
  XTX_inv = INVERT(XTX)
  XTY = XT # (signals[0:linindex, i])
  beta = XTX_inv # XTY
  pred = X # beta
  
  xsum = TOTAL((exposures[0:linindex] - avg(exposures[0:linindex]))^2)
  ysum = TOTAL((signals[0:linindex, i] - pred)^2)
  serr = sqrt(ysum/((n-2)*(xsum)))
  tstat = beta[1]/serr
  
  pval = 1 - t_pdf(tstat, n-2)
  
  ;----- Mark pixels defective based on nonlinearity and correlation (using a alpha level of 0.05) -----
  
  xc = i MOD xsize
  yc = FLOOR(i/xsize)
 
  inc = signals[linindex, i] - signals[0, i]
  
  ; If the pixel has a nonlinearity greater than 5%, a p-value larger than 0.05, or an weak linear slope...   
  if (nonlin gt 5) or (pval gt 0.05) or (inc lt (mx/2)) then begin
    mask[xc, yc] = 1
    print, '('+strtrim(string(xc),2)+', '+strtrim(string(yc),2)+')
  endif else begin
    mask[xc, yc] = 0
  endelse
  
  ; For troubleshooting/testing
  if nonlin gt 5 then begin
    nonlins[xc, yc] = nonlin
  endif
  if pval gt 0.05 then begin
    pvals[xc, yc] = pval
  endif
  if inc lt mx/2 then begin
    incs[xc, yc] = 1
  endif
endfor

end