; A program that tests the consistency of an illuminator and integrating sphere.
; Data is taken over a period of time at constant intervals and exposure, and the mean gray values are plotted against time.
; Ideally, a well-functioning illuminator should have a constant, zero-slope relationship. 

; This program uses the IDLAstro library.

dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

res = fltarr(n)
offset = 40

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

for i = 0, n-1 do begin
  img = readfits(list[i])
  res[i] = avg(img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset])
endfor

;stop           ; check for error images (after illuminator has shut off) 

t = FINDGEN(n, START=60, INCREMENT=0.25)
plot, t, res, psym=1, XTITLE='Time Elapsed (min)', YTITLE='Mean Signal [ADU]', $
  XMARGIN=[20,5], YMARGIN=[10,5], XRANGE=[min(t), max(t)], /YNOZERO

;======================== NONLINEARITY APPROACH ========================

; Plot a constant line that should represent perfect stability
YFIT = REPLICATE(avg(res), n)
OPLOT, t, YFIT, LINESTYLE=0

; Find the non-linearity (for the theoretical perfect stability)
mini = res[0]
maxi = 0.0
maxsig = 0.0
for index = 0, n-1 do begin
  curr = res[index] - avg(res)
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

;print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
print, 'Non-linearity (constant): ', nonlin*100, '%.'

; Plot the linear regression of the actual observed data (not perfect stability) and find the non-linearity
coeff = LINFIT(t[0:n-1], res[0:n-1])
YFIT = coeff[0] + coeff[1]*t
OPLOT, t, YFIT, LINESTYLE=1
print, 'Linear Regression: y =', coeff[1], 'x +', coeff[0]

mini0 = res[0]
maxi0 = 0.0
maxsig0 = 0.0
for index0 = 0, n-1 do begin
  curr0 = res[index0] - (coeff[1]*t[index0] + coeff[0])
  if curr0 gt maxi0 then begin
    maxi0 = curr0
  endif
  if curr0 lt mini0 then begin
    mini0 = curr0
  endif
  if res[index0] gt maxsig0 then begin
    maxsig0 = res[index0]
  endif
endfor

nonlin0 = (maxi0 + abs(mini0)) / maxsig0
;print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
print, 'Non-linearity (regression): ', nonlin0*100, '%.'

;; Plot a constant line that should represent perfect stability
;CFIT1 = REPLICATE(avg(res[0:129]), 130)
;OPLOT, t[0:129], CFIT1, LINESTYLE=0
;
;; Find the non-linearity (for the theoretical perfect stability)
;minic = res[0]
;maxic = 0.0
;maxcsig = 0.0
;for index = 0, 129 do begin
;  curr = res[index] - avg(res[0:129])
;  if curr gt maxic then begin
;    maxic = curr
;  endif
;  if curr lt minic then begin
;    minic = curr
;  endif
;  if res[index] gt maxcsig then begin
;    maxcsig = res[index]
;  endif
;endfor
;
;nonlinc = (maxic + abs(minic)) / maxcsig
;
;;print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
;print, 'Non-linearity ([0:129], constant): ', nonlinc*100, '%.'
;
;; Plot a constant line that should represent perfect stability
;CFIT2 = REPLICATE(avg(res[130:239]), 110)
;OPLOT, t[130:239], CFIT2, LINESTYLE=0
;
;; Find the non-linearity (for the theoretical perfect stability)
;minic1 = res[130]
;maxic1 = 0.0
;maxcsig1 = 0.0
;for index = 130, 239 do begin
;  curr = res[index] - avg(res[130:239])
;  if curr gt maxic1 then begin
;    maxic1 = curr
;  endif
;  if curr lt minic1 then begin
;    minic1 = curr
;  endif
;  if res[index] gt maxcsig1 then begin
;    maxcsig1 = res[index]
;  endif
;endfor
;
;nonlinc1 = (maxic1 + abs(minic1)) / maxcsig1
;
;;print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
;print, 'Non-linearity ([130:239], constant): ', nonlinc1*100, '%.'
;
;; Plot the linear regression of the actual observed data (not perfect stability) and find the non-linearity
;coeff1 = LINFIT(t[0:129], res[0:129])
;YFIT1 = coeff1[0] + coeff1[1]*t[0:129]
;OPLOT, t[0:129], YFIT1, LINESTYLE=1
;print, 'Linear Regression ([0:129]): y =', coeff1[1], 'x +', coeff1[0]
;
;mini0 = res[0]
;maxi0 = 0.0
;maxsig0 = 0.0
;for index0 = 0, 129 do begin
;  curr0 = res[index0] - (coeff1[1]*t[index0] + coeff1[0])
;  if curr0 gt maxi0 then begin
;    maxi0 = curr0
;  endif
;  if curr0 lt mini0 then begin
;    mini0 = curr0
;  endif
;  if res[index0] gt maxsig0 then begin
;    maxsig0 = res[index0]
;  endif
;endfor
;
;nonlin0 = (maxi0 + abs(mini0)) / maxsig0
;;print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
;print, 'Non-linearity ([0:129, regression): ', nonlin0*100, '%.'
;
;; Plot the linear regression of the actual observed data (not perfect stability) and find the non-linearity
;coeff2 = LINFIT(t[130:239], res[130:239])
;YFIT2 = coeff2[0] + coeff2[1]*t[130:239]
;OPLOT, t[130:239], YFIT2, LINESTYLE=1
;print, 'Linear Regression ([130:239]): y =', coeff2[1], 'x +', coeff2[0]
;
;mini1 = res[130]
;maxi1 = 0.0
;maxsig1 = 0.0
;for index1 = 130, 239 do begin
;  curr1 = res[index1] - (coeff2[1]*t[index1] + coeff2[0])
;  if curr1 gt maxi1 then begin
;    maxi1 = curr1
;  endif
;  if curr1 lt mini1 then begin
;    mini1 = curr1
;  endif
;  if res[index1] gt maxsig1 then begin
;    maxsig1 = res[index1]
;  endif
;endfor
;
;nonlin1 = (maxi1 + abs(mini1)) / maxsig1
;;print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
;print, 'Non-linearity ([130:239, regression): ', nonlin1*100, '%.'

end