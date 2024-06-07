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

num = FINDGEN(n)
plot, num, res, psym=1, XTITLE='Img Number', YTITLE='Mean Signal [ADU]', XMARGIN=[20,5], YMARGIN=[10,5], /YNOZERO

; Plot a constant line that should represent perfect stability
YFIT = REPLICATE(avg(res), n)
OPLOT, num, YFIT, LINESTYLE=0

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
coeff = LINFIT(num[0:n-1], res[0:n-1])
YFIT = coeff[0] + coeff[1]*num
OPLOT, num, YFIT, LINESTYLE=1
print, 'Linear Regression: y =', coeff[1], 'x +', coeff[0]

mini0 = res[0]
maxi0 = 0.0
maxsig0 = 0.0
for index0 = 0, n-1 do begin
  curr0 = res[index0] - (coeff[1]*num[index0] + coeff[0])
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

end

