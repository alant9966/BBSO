; /**
;  * niris_mask.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and identifies and outputs the coordinates of hot pixels and dead pixels.
;  * Coordinates of defect pixels are also stored in corresponding arrays for referencing.
;  * 
;  * This code functions the exact same as 'mask.pro', except accounting for scrambled data (images not in increasing exposure order).
;  * The images and darks must still be taken separately, such that the first half of the FITS files in the directory are the images and 
;  *    second half are the darks. If it is the other way around, simply switch the assignment of img and dk in the code.
;  *
;  * Only the image sorting algorithm will be commented on in this code. For detailed comments on the entire defect pixel identification
;  *    algorithm, refer to the 'mask.pro' code.
; */

dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

img0 = readfits(list[0], head)
xlen = head[3] & xsub = strmid(xlen, 0, strpos(xlen, '/')-2)
ylen = head[4] & ysub = strmid(ylen, 0, strpos(ylen, '/')-2)
xstart = strpos(xsub, ' ', /REVERSE_SEARCH) + 1
ystart = strpos(ysub, ' ', /REVERSE_SEARCH) + 1
xcount = strpos(xlen, '/') - xstart - 1
ycount = strpos(ylen, '/') - ystart - 1
xsize = float(strmid(xlen, xstart, xcount))
ysize = float(strmid(ylen, ystart, ycount))

exposures = fltarr(n/2)
signals = fltarr(n/2, (xsize*ysize)) & imgsig0 = signals
mask = fltarr(xsize, ysize)

nonlins = fltarr(xsize, ysize) & dead = nonlins & hot = nonlins
; To help in sorting the data, map each image and dark to its exposure, then signal value for every pixel
imgsig = fltarr(n/2, 1+(xsize*ysize)) & dksig = imgsig

;=============================== SORTING THE DATA ===============================

for i = 0, (n/2)-1 do begin
  ; Grab each image and its 'corresponding' dark, and their header files
  ;   (the dark will not actually correspond correctly, as the files are scrambled)
  dk = readfits(list[i+(n/2)], imghead)
  img = readfits(list[i], dkhead)

  ; Determine the exposures of the image and dark from the header files and add it to the corresponding array.
  ; Depending on the camera, the index for ...head[] may need to be changed, as the exposure may be listed in a slightly
  ;   different position in the header file. To do this, view the header file - each parameter listed is counted as one
  ;   index in the array. Count from the top until the index of the exposure parameter is determined (start from index 0).
  exp1 = imghead[24] & expsub1 = strmid(exp1, 0, strpos(exp1, '.'))
  start1 = strpos(expsub1, ' ', /REVERSE_SEARCH) + 1
  count1 = strpos(exp1, '/') - start1 - 1
  imgsig[i,0] = float(strmid(exp1, start1, count1))

  exp2 = dkhead[24] & expsub2 = strmid(exp2, 0, strpos(exp2, '.'))
  start2 = strpos(expsub2, ' ', /REVERSE_SEARCH) + 1
  count2 = strpos(exp2, '/') - start2 - 1
  dksig[i,0] = float(strmid(exp2, start2, count2))

  ; Choose the middle frame of the image and dark...
  img0 = img[*,*, FLOOR((size(img))[3]/2)]
  dk0 = dk[*,*, FLOOR((size(dk))[3]/2)]
 
  ; ...then add the signals for each pixel to the corresponding array.
  curr = 1UL;
  for y = 0, ysize-1 do begin
    for x = 0, xsize-1 do begin
      imgsig[i,curr] = img0[x,y]
      dksig[i,curr] = dk0[x,y]
      curr++
    endfor
  endfor
endfor

; Sort the images and darks in order of increasing exposure 
;   (the resulting array lists the indexes of the exposures in increasing order).
imgsort = sort(imgsig[*,0])
dksort = sort(dksig[*,0])

for i = 0, (n/2)-1 do begin
  ; The indexes of the image and dark with the same exposure
  imgi = imgsort[i]
  dki = dksort[i]

  ; Subtract the dark (for each pixel) and add the resulting value to the signals array, which is now sorted
  for j = 0, (xsize*ysize)-1 do begin
    if (FLOAT(imgsig[imgi,j+1]) - FLOAT(dksig[dki,j+1])) lt 0 then begin
      signals[i,j] = 0
    endif else begin
      signals[i,j] = FLOAT(imgsig[imgi,j+1]) - FLOAT(dksig[dki,j+1])
    endelse
    imgsig0[i,j] = imgsig[imgi,j+1]
  endfor
  
  ; Add the corresponding exposure to the exposures array
  exposures[i] = imgsig[imgi,0]
endfor

;================================================================================

READ, linrange, PROMPT='Please enter the linear range (ADU): '
linindex = (n/2)-1
for i = 0, (n/2)-2 do begin
  avgsig = 0.0
  avgsignext = 0.0
  count = 0.0
  for k = 0, xsize*ysize-1 do begin
    if ~(avg(signals[*,k]) lt 10000) then begin
      avgsig += signals[i,k]
      avgsignext += signals[i+1,k]
      count++
    endif
  endfor
  avgsig /= count
  avgsignext /= count

  if (avgsig le linrange) and (avgsignext ge linrange) then begin
    linindex = i
    break
  endif
  print, i
endfor

interval = exposures[1] - exposures[0]

print, 'Coordinates of nonlinear pixels: '

for i = 0, (xsize*ysize)-1 do begin
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

  if maxsig eq 0 then begin
    nonlin = 100
  endif else begin
    nonlin = ((maxi+abs(mini))/maxsig)*100
  endelse

  xc = i MOD xsize
  yc = FLOOR(i/xsize)

  if (nonlin gt 15) then begin
    mask[xc, yc] = 1
    nonlins[xc, yc] = nonlin
    print, '('+strtrim(string(xc),2)+', '+strtrim(string(yc),2)+')
  endif else begin
    mask[xc, yc] = 0
  endelse 
endfor

print, 'Coordinates of dead pixels and hot pixels: '

hotindex = FLOOR(linindex/2)
for i = 0, (xsize*ysize)-1 do begin
  xc = i MOD xsize
  yc = FLOOR(i/xsize)

  if imgsig0[(n/2)-1,i] lt 500 then begin
    dead[xc,yc] = 1
    mask[xc,yc] = 1
    print, 'dead: ('+strtrim(string(xc),2)+', '+strtrim(string(yc),2)+')
  endif

  if imgsig0[hotindex,i] ge avg(signals[(n/2)-1,i]) then begin
    hot[xc,yc] = 1
    mask[xc,yc] = 1
    print, 'hot: ('+strtrim(string(xc),2)+', '+strtrim(string(yc),2)+')
  endif
endfor

end
