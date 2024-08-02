; /**
;  * niris_regression.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and plots the exposure (ms) vs. mean signal (ADU) curve.
;  * Linear regression is performed over the desired linear range, and nonlinearity is calculated.
;  * Refer to the console for detailed information about the regression line and nonlinearity value.
;  *
;  * This code functions the exact same as 'nonlinearity.pro', except accounting for scrambled data (images not in increasing exposure order)
;  *    and accounting for the use of a non-uniform light source such as the sun for testing.
;  * The images and darks must still be taken separately, such that the first half of the FITS files in the directory are the images and
;  *    second half are the darks. If it is the other way around, simply switch the assignment of img and dk in the code.
;  *
;  * For detailed comments on the entire linear regression algorithm, refer to the 'nonlinearity.pro' code. Only major differences in the
;  *    code will be commented on in this file.
; */

dir = dialog_pickfile(/directory,path='c:\')
list = file_search(dir + '*.fts' , count = n)

offset = 20

exposures = fltarr(n/2) & res = exposures
signals = fltarr(n/2, (2*offset+1)*(2*offset+1))
imgsig = fltarr(n/2, 1+((2*offset+1)*(2*offset+1))) & dksig = imgsig

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

;=============================== SORTING THE DATA ===============================

for i = 0, (n/2)-1 do begin
  ; Grab each image and its 'corresponding' dark, and their header files
  ;   (the dark will not actually correspond correctly, as the files are scrambled)
  img = readfits(list[i+(n/2)], imghead)
  dk = readfits(list[i], dkhead)
  
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
  
  ; Choose the middle frame of the image and dark and identify the ROI...
  img0 = img[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, FLOOR((size(img))[3]/2)]
  dk0 = dk[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, FLOOR((size(dk))[3]/2)]
  
  ; ...then add the signals for each pixel to the corresponding array.
  curr = 1UL;
  for y = 0, (2*offset) do begin
    for x = 0, (2*offset) do begin
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
  for j = 0, (2*offset+1)*(2*offset+1)-1 do begin
    if (FLOAT(imgsig[imgi,j+1]) - FLOAT(dksig[dki,j+1])) lt 0 then begin
      signals[i,j] = 0
    endif else begin
      signals[i,j] = FLOAT(imgsig[imgi,j+1]) - FLOAT(dksig[dki,j+1])
    endelse
  endfor
  
  ; Add the corresponding exposure to the exposures array
  exposures[i] = imgsig[imgi,0]
endfor

; Add the sorted mean exposures to the res array
for i = 0, (n/2)-1 do begin
  res[i] = avg(signals[i,*])
endfor

;================================================================================

plot, exposures, res, psym=1, XTITLE='Exposure Time [ms]', YTITLE='Mean Signal [ADU]'

cont = 'yes'
while cont eq 'yes' do begin
  READ, cilms, PROMPT='Please input the linear max (ADU): '

  ind = 0
  for i = 0, (n/2)-2 do begin
    if (res[i] le cilms) and (res[i+1] ge cilms) then begin
      ind = i
      break
    endif
  endfor

  coeff = LINFIT(exposures[0:ind], res[0:ind])
  YFIT = coeff[0] + coeff[1]*exposures

  OPLOT, exposures, YFIT
  print, 'Linear Regression: y = '+strtrim(string(coeff[1]),2)+'x + '+strtrim(string(coeff[0]),2)

  mini = res[ind]
  maxi = 0.0
  maxsig = 0.0
  for index = 0, ind do begin
    curr = res[index] - (coeff[1]*exposures[index] + coeff[0])
    if curr gt maxi then begin
      maxi = curr
    endif
    if curr lt mini then begin
      mini = curr
    endif
    if res[index] gt maxsig then begin
      maxsig = res[index]
    endif
  endfor

  nonlin = (maxi + abs(mini)) / maxsig

  print, 'The maximum deviation is '+strtrim(string(maxi),2)+', the minimum deviation is '+strtrim(string(mini),2)+' and the maximum signal is '+strtrim(string(maxsig),2)+'.'
  print, 'The nonlinearity over the linear range is '+strtrim(string(nonlin*100),2)+'%.'

  str = ''
  READ, str, PROMPT='Continue regression? (y/n) '
  if (str eq 'n') then begin
    cont = 'no'
  endif
endwhile

end