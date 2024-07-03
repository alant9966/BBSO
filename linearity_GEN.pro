; Prompts the user to select a folder, and reads all FITS files into an array of file names
dir = dialog_pickfile(/directory,path='c:\')
list = file_search(dir + '*.fts' , count = n)

;; Prompts the user to input the data's starting exposure and interval
;READ, t0, PROMPT='Please input the starting exposure (in ms): '
;READ, interval, PROMPT='Please input the exposure interval: '

exposures = fltarr(n/2)                
signals = fltarr(n/2) & darks = signals 

;offset = 20
offset = 40

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

; If the dark field images are grouped at the end...
for i = 0, (n/2)-1 do begin
  img = readfits(list[i], header)
  dk = readfits(list[i+(n/2)])
  
  signals[i] = avg(img[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, FLOOR((size(img))[3]/2)] )
  darks[i] = avg(dk[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, FLOOR((size(dk))[3]/2)] )
  
  ; The recorded exposure for each image is added to the exposures array.
  exp = header[26] & expsub = strmid(exp, 0, strpos(exp, '.'))
  start = strpos(expsub, ' ', /REVERSE_SEARCH) + 1
  count = strpos(exp, '/') - start - 1
  exposures[i] = float(strmid(exp, start, count))
endfor

; If the images and dark fields are taken together...
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

res = signals - darks
interval = exposures[1] - exposures[0]

;stop

; Plot the data and perform a linear regression
plot, exposures, res, psym=1, XTITLE='Exposure Time [ms]', YTITLE='Mean Signal [ADU]'

; Allows the user to test linear regression intervals until satisfied
;cont = boolean(1)
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
  print, 'Equation of linear regression from 0 to '+strtrim(string(cilms),2)+' ADU: y = '+strtrim(string(coeff[1]),2)+'x + '+strtrim(string(coeff[0]),2)
  
  ; Determine the nonlinearity of the data
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