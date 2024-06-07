dir = dialog_pickfile(/directory,path='e:\')
list = file_search( dir + '*.fts' , count = n )

r1 = fltarr(100) & d1 = r1
r2 = fltarr(100) & d2 = r2
rvar = fltarr(100) & dvar = rvar

; plot var(sig) - var(dark) vs mean(sig) - mean(dark) instead
 
offset = 40                                   

for i = 0, 99 do begin
  img = readfits(list[i])
  r1[i] = avg( img[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 0] )                   
  r2[i] = avg( img[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 4] )                   
  
  img1 = img[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 0]
  img2 = img[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 4]
  rvar[i] = 0.5*(avg((img1-img2)^2))
endfor

for j = 100, 109 do begin
  dark = readfits(list[j])
  for k = 0, 9 do begin
    d1[10*(j-100)+k] = avg( dark[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 0] )
    d2[10*(j-100)+k] = avg( dark[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 4] )
    
    dimg1 = dark[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 0]
    dimg2 = dark[1024-256-offset : 1024-256+offset, 1024+256-offset : 1024+256+offset, 4]
    dvar[k] = 0.5*(avg((dimg1-dimg2)^2))
  endfor
endfor

r1 = r1 - d1 & r2 = r2 - d2
sig = (r1+r2)/2

var = rvar - dvar
coeff = LINFIT(sig[0:67], var[0:67])  ; 1.7 ms
YFIT = coeff[0] + coeff[1]*sig
plot, sig, var, psym=1, XTITLE='Mean Signal [ADU]', YTITLE='Temporal Variance [ADU]'

;OPLOT, sig, YFIT
print, 'Equation of the linear fit: y = ', coeff[1], 'x + ', coeff[0]

end