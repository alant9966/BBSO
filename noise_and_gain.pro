; /**
;  * noise_and_gain.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and plots the photon-transfer curve (mean signal (ADU) vs. signal variance (ADU)).
;  * Linear regression is performed over the desired linear range, and an estimated readout noise and gain is calculated from the resulting equation.
;  * Refer to the console for detailed information about the readout noise and gain values.
;  *
;  * This code assumes that:
;  *    1. The darks and images are taken separately, and each in order of increasing exposure
;  *        i.e. an example directory would be (1.5 ms img, 3 ms img, 4.5 ms img, 1.5 ms dark, 3 ms dark, 4.5 ms dark) from top to bottom
;  *    2. A uniform light source is used, such that every pixel on the camera is equally and constantly illuminated.
;  *
;  * If one or both of these criteria are not met (if the sun is used instead of uniform illumination), use the "niris_noise.pro" code instead,
;  * which ignores irrelevant pixels (outside of the focus area) and accounts for scrambled data (images not in increasing exposure order).
; */

; Select the directory with the FITS image files
dir = dialog_pickfile(/directory,path='e:\')
list = file_search( dir + '*.fts' , count = n )

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

; 40 x 40 pixel ROI
offset = 20

; Create arrays for the images, darks, and variances
d1 = fltarr((2*offset)+1,(2*offset)+1) & r1 = d1
d2 = d1 & r2 = d1 & z = d1
dev = fltarr(n/2) & var = dev
sig = fltarr(n/2)

; For each image...
for i = 0, n/2 - 1 do begin
  ; Grab the image and its corresponding dark
  img = readfits(list[i])
  dk = readfits(list[i+(n/2)])

  ; Subtract darks for the first and last frames of the image
  r1 = float(img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 0])
  d1 = float(dk[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 0])
  r1 = r1 - d1
  
  r2 = float(img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, (size(img))[3]-1])
  d2 = float(dk[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, (size(img))[3]-1])
  r2 = r2 - d2
  
  ; *Not exactly sure what zavg and zdev represent but they are necessary for final calculations
  z = r2/r1
  zavg = avg(z)
  zdev = stddev(z)

  ; Average the two frames
  r0 = (r1+r2)/2.                           
  sig[i] = avg(r0) ;& print,s[i]                

  ; Calculate the standard deviation and variance, adding each to their corresponding arrays.
  ; Note that the standard deviation is not used unless plotting the alternative version of the photon-transfer curve (log signal vs. log std dev).
  dev[i] = sig[i]/zavg * zdev/sqrt(2.)                 
  var[i] = (sig[i]/zavg * zdev/sqrt(2.))^2              
endfor

; Plot the photon transfer curve (mean signal vs. signal variance)
plot, sig, var, psym=1, XTITLE='Mean Signal [ADU]', YTITLE='Temporal Variance [ADU]', $
  XMARGIN=[20,5], YMARGIN=[8,3];, MAX_VALUE=1000000.0;, XRANGE=[0,70000], YRANGE=[0,15000]
  
; Alternatively, you can plot the (log signal vs. log standard deviation) curve, but this may require slight different calculations for gain and noise.
;plot, sig, dev, psym=1, XTITLE='Log Mean Signal [ADU]', YTITLE='Log Noise [ADU]', $
;  XMARGIN=[20,5], YMARGIN=[8,3], MAX_VALUE=1000.0, /XLOG, /YLOG;, XRANGE=[0,70000], YRANGE=[0,15000]

; Determines the index to regress over (the linear range)
READ, linr, PROMPT='Please input the linear max (ADU): '
ind = 0
for i = 0, (n/2)-2 do begin
  if (sig[i] le linr) and (sig[i+1] ge linr) then begin
    ind = i
    break
  endif
endfor

; Plot the linear regression and output the equation the console
coeff = LINFIT(sig[0:ind], var[0:ind])
YFIT = coeff[0] + coeff[1]*sig

OPLOT, sig, YFIT
print, 'Linear Regression: y = '+strtrim(string(coeff[1]),2)+'x + '+strtrim(string(coeff[0]),2)

; Calculate the gain and noise from the slope and intercept of the regression
gain = 1/coeff[1]
rnoise = sqrt(coeff[0]*(gain^2))
ernoise = rnoise * gain

print, 'Gain: '+strtrim(string(gain),2)+' e-/ADU'
print, 'Readout Noise: '+strtrim(string(rnoise),2)+' ADU ('+strtrim(string(ernoise),2)+' e-)'

end