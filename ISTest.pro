; TODO: Test the performance of the integrating sphere - run the sphere for around 10 minutes and graph the 
; change in light intensity over time to make sure that it is constant and stable
; 
; Fixed exposure / illuminator intensity, take an image every x seconds and plot the signal in ADU
; Can also test the uniformity by graphing the entire lens, then comparing graphs of several portions
; for uniformity in intensity and stability
 
dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

res = fltarr(n)
offset = 20

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

for i = 0, n-1 do begin
  img = readfits(list[i])
  res[i] = avg(img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset])
endfor

num = FINDGEN(n)
plot, num, res, psym=1, XTITLE='Img Number', YTITLE='Mean Signal [ADU]'

; Plot a constant line that should represent perfect stability
YFIT = avg(res)
OPLOT, num, YFIT, LINESTYLE=0

; Find the non-linearity
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

coeff = LINFIT(num[0:n-1], res[0:n-1])
YFIT = coeff[0] + coeff[1]*num

OPLOT, exposures, YFIT, LINESTYLE=1
print, 'Linear Regression: y =', coeff[1], 'x +', coeff[0]

end

