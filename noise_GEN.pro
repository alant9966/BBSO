dir = dialog_pickfile(/directory,path='e:\')
list = file_search( dir + '*.fts' , count = n )

r1 = fltarr(n/2) & d1 = r1
r2 = fltarr(n/2) & d2 = r2
rvar = fltarr(n/2) & dvar = rvar

offset = 20

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

; If the dark field images are grouped at the end...
for i = 0, (n/2)-1 do begin
  img = readfits(list[i])
  
  img1 = img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 0]
  img2 = img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 4]
  
  r1[i] = avg( img1 )
  r2[i] = avg( img2 )
  
  rvar[i] = 0.5*(avg((img1-img2)^2))
endfor

for j = n/2, n-1 do begin
  dark = readfits(list[j])
  
  dimg1 = dark[xcoord-offset : xcoord+offset, xcoord-offset : xcoord+offset, 0]
  dimg2 = dark[xcoord-offset : xcoord+offset, xcoord-offset : xcoord+offset, 4]
  
  d1[j-(n/2)] = avg( dimg1 )
  d2[j-(n/2)] = avg( dimg2 )

  dvar[j-(n/2)] = 0.5*(avg((dimg1-dimg2)^2))
endfor

; If the images and dark fields are taken together...
;for i = 0, n/2 - 1 do begin
;  img = readfits(list[2*i+1])
;  dark = readfits(list[2*i])
;  
;  img1 = img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 0]
;  img2 = img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 4]
;  dimg1 = dark[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 0]
;  dimg2 = dark[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 4]
;
;  r1[i] = avg( img1 )
;  r2[i] = avg( img2 )
;  d1[i] = avg( dimg1 )
;  d2[i] = avg( dimg2 )
;
;  rvar[i] = 0.5*(avg((img1-img2)^2))
;  dvar[i] = 0.5*(avg((dimg1-dimg2)^2))
;endfor

r1 = r1 - d1 & r2 = r2 - d2
sig = (r1+r2)/2

var = rvar - dvar
plot, sig, var, psym=1, XTITLE='Mean Signal [ADU]', YTITLE='Temporal Variance [ADU]', $
  XRANGE=[0,70000], YRANGE=[0,15000], XMARGIN=[20,5], YMARGIN=[8,3]

READ, cil, PROMPT='Please input the linear max (70% saturation, in ADU): '

; will give an out-of-bounds error if the loop reaches the end of the array (and calls k+1),
; but this should never happen as the input should be 70% of saturation anyway
ind = 0
for k = 0, n/2 do begin
  if (sig[k] eq cil) or (sig[k] lt cil and sig[k+1] gt cil) then begin
    ind = k
    break
  endif
endfor

coeff = LINFIT(sig[0:ind], var[0:ind])
YFIT = coeff[0] + coeff[1]*sig

OPLOT, sig, YFIT
print, 'Linear Regression: y =', coeff[1], 'x +', coeff[0]

end