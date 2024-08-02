; /**
;  * noise_and_gain.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and plots the photon-transfer curve (mean signal (ADU) vs. signal variance (ADU)).
;  * Linear regression is performed over the desired linear range, and an estimated readout noise and gain is calculated from the resulting equation.
;  * Refer to the console for detailed information about the readout noise and gain values.
;  *
;  * This code functions the exact same as 'noise_and_gain.pro', except accounting for scrambled data (images not in increasing exposure order)
;  *    and accounting for the use of a non-uniform light source such as the sun for testing.
;  * The images and darks must still be taken separately, such that the first half of the FITS files in the directory are the images and
;  *    second half are the darks. If it is the other way around, simply switch the assignment of img and dk in the code.
;  *
;  * The method/algorithm used to sort the scrambled data is the exact same as for the other 'niris____.pro' scripts.
; */

dir = dialog_pickfile(/directory,path='e:\')
list = file_search( dir + '*.fts' , count = n )

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

offset = 20

d1 = fltarr((2*offset)+1,(2*offset)+1) & r1 = d1
d2 = d1 & r2 = d1 & z = d1
dev = fltarr(n/2) & var = dev
sig = fltarr(n/2)

signals1 = fltarr(n/2, (2*offset+1)*(2*offset+1)) & signals2 = signals1 & signals = signals1
imgsig1 = fltarr(n/2, 1+((2*offset+1)*(2*offset+1))) & dksig1 = imgsig1
imgsig2 = imgsig1 & dksig2 = imgsig1

for i = 0, (n/2)-1 do begin
  img = readfits(list[i+(n/2)], imghead)
  dk = readfits(list[i], dkhead)

  exp1 = imghead[24] & expsub1 = strmid(exp1, 0, strpos(exp1, '.'))
  start1 = strpos(expsub1, ' ', /REVERSE_SEARCH) + 1
  count1 = strpos(exp1, '/') - start1 - 1
  imgsig[i,0] = float(strmid(exp1, start1, count1))

  exp2 = dkhead[24] & expsub2 = strmid(exp2, 0, strpos(exp2, '.'))
  start2 = strpos(expsub2, ' ', /REVERSE_SEARCH) + 1
  count2 = strpos(exp2, '/') - start2 - 1
  dksig[i,0] = float(strmid(exp2, start2, count2))

  r1 = img[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, 0]
  r2 = img[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, (size(img))[3]-1]
  d1 = dk[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, 0]
  d2 = dk[xcoord-offset:xcoord+offset, ycoord-offset:ycoord+offset, (size(dk))[3]-1]

  curr = 1UL;
  for y = 0, (2*offset) do begin
    for x = 0, (2*offset) do begin
      imgsig1[i,curr] = r1[x,y]
      imgsig2[i,curr] = r2[x,y]
      dksig1[i,curr] = d1[x,y]
      dksig2[i,curr] = d2[x,y]
      curr++
    endfor
  endfor
endfor

;stop

imgsort = sort(imgsig[*,0])
dksort = sort(dksig[*,0])

for i = 0, (n/2)-1 do begin
  imgi = imgsort[i]
  dki = dksort[i]

  for j = 0, (2*offset+1)*(2*offset+1)-1 do begin
    if (FLOAT(imgsig1[imgi,j+1]) - FLOAT(dksig1[dki,j+1])) lt 0 then begin
      signals1[i,j] = 0
    endif else begin
      signals1[i,j] = FLOAT(imgsig1[imgi,j+1]) - FLOAT(dksig1[dki,j+1])
    endelse
    
    if (FLOAT(imgsig2[imgi,j+1]) - FLOAT(dksig2[dki,j+1])) lt 0 then begin
      signals2[i,j] = 0
    endif else begin
      signals2[i,j] = FLOAT(imgsig2[imgi,j+1]) - FLOAT(dksig2[dki,j+1])
    endelse
  endfor
endfor

stop

for i = 0, (n/2)-1 do begin
  for j = 0, (2*offset+1)*(2*offset+1)-1 do begin
    signals[i,j] = (signals1[i,j]+signals2[i,j])/2.
  endfor
endfor

for i = 0, n/2 - 1 do begin
  r01 = signals1[i,*]
  r02 = signals2[i,*]

  z = r02/r01
  zavg = avg(z)
  zdev = stddev(z)

  sig[i] = avg(signals[i,*]) ;& print,s[i]

  dev[i] = sig[i]/zavg * zdev/sqrt(2.)
  var[i] = (sig[i]/zavg * zdev/sqrt(2.))^2
  ;print, i
endfor

plot, sig, var, psym=1, XTITLE='Mean Signal [ADU]', YTITLE='Temporal Variance [ADU]', $
  XMARGIN=[20,5], YMARGIN=[8,3];, MAX_VALUE=1000000.0;, XRANGE=[0,70000], YRANGE=[0,15000]
;plot, sig, dev, psym=1, XTITLE='Log Mean Signal [ADU]', YTITLE='Log Noise [ADU]', $
;  XMARGIN=[20,5], YMARGIN=[8,3], MAX_VALUE=1000.0, /XLOG, /YLOG;, XRANGE=[0,70000], YRANGE=[0,15000]

READ, linr, PROMPT='Please input the linear max (ADU): '
ind = 0
for i = 0, (n/2)-2 do begin
  if (sig[i] le linr) and (sig[i+1] ge linr) then begin
    ind = i
    break
  endif
endfor
coeff = LINFIT(sig[0:ind], var[0:ind])
YFIT = coeff[0] + coeff[1]*sig

OPLOT, sig, YFIT
print, 'Linear Regression: y = '+strtrim(string(coeff[1]),2)+'x + '+strtrim(string(coeff[0]),2)

gain = 1/coeff[1]
rnoise = sqrt(coeff[0]*(gain^2))
ernoise = rnoise * gain

print, 'Gain: '+strtrim(string(gain),2)+' e-/ADU'
print, 'Readout Noise: '+strtrim(string(rnoise),2)+' ADU ('+strtrim(string(ernoise),2)+' e-)'

end
