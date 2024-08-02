; /**
;  * nonlinearity.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and plots the exposure (ms) vs. mean signal (ADU) curve.
;  * Linear regression is performed over the desired linear range, and nonlinearity is calculated.
;  * Refer to the console for detailed information about the regression line and nonlinearity value.
;  *
;  * This code assumes that:
;  *    1. The darks and images are taken separately, and each in order of increasing exposure
;  *        i.e. an example directory would be (1.5 ms img, 3 ms img, 4.5 ms img, 1.5 ms dark, 3 ms dark, 4.5 ms dark) from top to bottom
;  *    2. A uniform light source is used, such that every pixel on the camera is equally and constantly illuminated.
;  *
;  * If one or both of these criteria are not met (if the sun is used instead of uniform illumination), use the "niris_regression.pro" code instead,
;  * which ignores irrelevant pixels (outside of the focus area) and accounts for scrambled data (images not in increasing exposure order).
; */

; Select the directory with the image FITS files
dir = dialog_pickfile(/directory,path='c:\')
list = file_search(dir + '*.fts' , count = n)

exposures = fltarr(n/2)                
signals = fltarr(n/2) & darks = signals 

; Sets the ROI to a 40 x 40 pixel area. Can be freely changed to any desired ROI dimension value.
offset = 20

; e.g. if the desired ROI is the center of the image and the aperture dimensions are 2048 x 2048, input
;   1024 for both prompts (as (1024, 1024) is the center coordinate).
print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

; If the images and darks are taken separately...
for i = 0, (n/2)-1 do begin  
  ; Grab each image and its corresponding dark (same exposures)
  img = readfits(list[i], header)
  idk = readfits(list[i+(n/2)])
  
  ; Add the average signal of each's ROI to their corresponding arrays
  signals[i] = avg(img[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, FLOOR((size(img))[3]/2)] )
  darks[i] = avg(dk[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, FLOOR((size(dk))[3]/2)] )
  
  ; The recorded exposure for each image is added to the exposures array
  exp = header[26] & expsub = strmid(exp, 0, strpos(exp, '.'))
  start = strpos(expsub, ' ', /REVERSE_SEARCH) + 1
  count = strpos(exp, '/') - start - 1
  exposures[i] = float(strmid(exp, start, count))
endfor

; If the images and darks are taken together/alternating...
;for i = 0, n/2 - 1 do begin
;  img = readfits(list[2*i+1])
;  dk = readfits(list[2*i])
;  signals[i] = avg( img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, FLOOR((size(img))[3]/2)] )
;  darks[i] = avg( dk[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, FLOOR((size(img))[3]/2)] )
;  
;  xp = header[26] & expsub = strmid(exp, 0, strpos(exp, '.'))
;  start = strpos(expsub, ' ', /REVERSE_SEARCH) + 1
;  count = strpos(exp, '/') - start - 1
;  exposures[i] = float(strmid(exp, start, count))
;endfor

; Subtract the darks from the images
res = signals - darks
interval = exposures[1] - exposures[0]

;stop

; Plot the data (exposure vs. mean signal0
plot, exposures, res, psym=1, XTITLE='Exposure Time [ms]', YTITLE='Mean Signal [ADU]'

; Allows the user to test linear regression intervals until satisfied
cont = 'yes'
while cont eq 'yes' do begin
  READ, cilms, PROMPT='Please input the linear max (ADU): '
  
  ; Determines the index of the signal value closest to the linear max
  ind = 0
  for i = 0, (n/2)-2 do begin
    if (res[i] le cilms) and (res[i+1] ge cilms) then begin
      ind = i
      break
    endif
  endfor
  
  ; Calculate and plot the linear regression
  coeff = LINFIT(exposures[0:ind], res[0:ind])
  YFIT = coeff[0] + coeff[1]*exposures
  
  OPLOT, exposures, YFIT
  print, 'Equation of linear regression from 0 to '+strtrim(string(cilms),2)+' ADU: y = '+strtrim(string(coeff[1]),2)+'x + '+strtrim(string(coeff[0]),2)
  
  ; Calculate the nonlinearity of the data
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