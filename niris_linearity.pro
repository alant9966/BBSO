dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

t0 = 0.025
interval = 0.025

t = ( FINDGEN(100) ) * interval + t0          ; [0.025, 0.05, 0.075, ... , 2.5] X-AXIS
d = fltarr(100) & r = d

offset = 40                                   ; select a 40 x 40 area to examine

for i = 0, 99 do begin                                                                    ; computes the average signal value of the center 40 x 40
  img = readfits(list[i])                                                                 ;     portion of the array for frame 3
  r[i] = avg( img[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 3] )              ; *would it be better to compute the avg of all five frames?
endfor

for j = 100, 109 do begin                                                                 ; does the same, but for the dark fields
  dark = readfits(list[j])
  for k = 0, 9 do begin                                                                         ; since only 10 dark fields were taken, every 10 images
    d[10*(j-100)+k] = avg( dark[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 3] )      ; calibrate with the same dark field
   endfor
endfor

res = r - d                                     ; final result array of the mean signal values

; plot the data and perform a linear fit

; DARK SUBTRACTION
coeff = LINFIT(t[0:59], res[0:59], CHISQR=chi)  ; assumes that the linear range is from 0.25 ms to 1.7 ms
YFIT = coeff[0] + coeff[1]*t
plot, t, res, psym=1, XTITLE='Exposure Time [ms]', YTITLE='Mean Signal [ADU]'       

OPLOT, t, YFIT
print, 'Equation of the linear regression from x = 0.025 ms to x = 1.7 ms: y = ', coeff[1], 'x + ', coeff[0]

;stop

; find the minimum and maximum deviations of the data (from the linear regression line)
mini = res[59]
maxi = 0.0
maxsig = 0.0
for index = 0, 59 do begin
  curr = res[index] - (coeff[1]*t[index] + coeff[0])
  if curr gt maxi then begin
    maxi = curr 
  endif
  if curr lt mini then begin
    mini = curr
  endif
  if res[index] gt maxsig then begin
    maxsig = res[index]
  endif  
  ;print, 'curr: ', curr, ' min: ', mini, ', max: ', maxi, ', maxsig: ', maxsig
endfor

; NO DARK SUBTRACTION
;coeff = LINFIT(t[0:71], r[0:71], CHISQR=chi)  ; assumes that the linear range is from 0.25 ms to 1.7 ms
;YFIT = coeff[0] + coeff[1]*t
;plot, t, r, psym=1, XTITLE='Exposure Time [ms]', YTITLE='Mean Signal [ADU]'
;OPLOT, t, YFIT
;
;print, 'Equation of the linear regression from x = 0.025 ms to x = 1.8 ms: y = ', coeff[1], 'x + ', coeff[0]
;
;; find the minimum and maximum deviations of the data (from the linear regression line)
;mini = r[72]
;maxi = 0.0
;maxsig = 0.0
;for index = 0, 71 do begin
;  curr = r[index] - (1177.24*t[index] + 124.645)
;  if curr gt maxi then begin
;    maxi = curr
;  endif
;  if curr lt mini then begin
;    mini = curr
;  endif
;  if r[index] gt maxsig then begin
;    maxsig = r[index]
;  endif
;  ;print, 'curr: ', curr, ' min: ', mini, ', max: ', maxi, ', maxsig: ', maxsig
;endfor


nonlin = (maxi + abs(mini)) / maxsig

print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
print, 'The nonlinearity is ', nonlin*100, '%.'

end