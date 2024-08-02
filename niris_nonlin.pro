; /**
;  * niris_nonlin.pro
;  * Alan Tong
;  *
;  * Reads a directory of FITS files and plots nonlinearity distribution histograms for all pixels at several potential linear ranges.
;  * Histograms will be saved as .jpg image files in a chosen directory. The console will read "Finished." when the process is finished.
;  * The user may use the histograms to determine what the ideal linear range is to create the nonlinearity curve.
;  * Then, run the 'niris_regression.pro' file on the data to plot the nonlinearity curve.
;  *
;  * This code functions the exact same as 'linear_range.pro', except accounting for scrambled data (images not in increasing exposure order)
;  *    and accounting for the use of a non-uniform light source such as the sun for testing.
;  * The images and darks must still be taken separately, such that the first half of the FITS files in the directory are the images and
;  *    second half are the darks. If it is the other way around, simply switch the assignment of img and dk in the code.
;  *
;  * For detailed comments on the entire linear regression algorithm, refer to the 'linear_range.pro' code. Only major differences in the
;  *    code will be commented on in this file.
; */

dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

print, 'Initializing parameters and variables...'
img0 = readfits(list[0], head)
xlen = head[3] & xsub = strmid(xlen, 0, strpos(xlen, '/')-2)
ylen = head[4] & ysub = strmid(ylen, 0, strpos(ylen, '/')-2)
xstart = strpos(xsub, ' ', /REVERSE_SEARCH) + 1
ystart = strpos(ysub, ' ', /REVERSE_SEARCH) + 1
xcount = strpos(xlen, '/') - xstart - 1
ycount = strpos(ylen, '/') - ystart - 1
xsize = float(strmid(xlen, xstart, xcount))
ysize = float(strmid(ylen, ystart, ycount))

exposures = fltarr(n/2)
signals = fltarr(n/2, (xsize*ysize))
nonlins = fltarr(xsize*ysize)

; To help in sorting the data, map each image and dark to its exposure, then signal value for every pixel
imgsig = fltarr(n/2, 1+(xsize*ysize)) & dksig = imgsig

;=============================== SORTING THE DATA ===============================

for i = 0, (n/2)-1 do begin
  ; Grab each image and its 'corresponding' dark, and their header files
  ;   (the dark will not actually correspond correctly, as the files are scrambled)
  dk = readfits(list[i+(n/2)], imghead)
  img = readfits(list[i], dkhead)
  
  ; Determine the exposures of the image and dark from the header files and add it to the corresponding array.
  ; Depending on the camera, the index for ...head[] may need to be changed, as the exposure may be listed in a slightly
  ;   different position in the header file. To do this, view the header file - each parameter listed is counted as one
  ;   index in the array. Count from the top until the index of the exposure parameter is determined (start from index 0).
  exp1 = imghead[24] & expsub1 = strmid(exp1, 0, strpos(exp1, '.'))
  start1 = strpos(expsub1, ' ', /REVERSE_SEARCH) + 1
  count1 = strpos(exp1, '/') - start1 - 1
  imgsig[i,0] = float(strmid(exp1, start1, count1))

  exp2 = dkhead[24] & expsub2 = strmid(exp2, 0, strpos(exp2, '.'))
  start2 = strpos(expsub2, ' ', /REVERSE_SEARCH) + 1
  count2 = strpos(exp2, '/') - start2 - 1
  dksig[i,0] = float(strmid(exp2, start2, count2))
  
  ; Choose the middle frame of the image and dark...
  img0 = img[*,*, FLOOR((size(img))[3]/2)]
  dk0 = dk[*,*, FLOOR((size(dk))[3]/2)]
  
  ; ...then add the signals for each pixel to the corresponding array.
  curr = 1UL;
  for y = 0, ysize-1 do begin
    for x = 0, xsize-1 do begin
      imgsig[i,curr] = img0[x,y]
      dksig[i,curr] = dk0[x,y]
      curr++
    endfor
  endfor
endfor

; Sort the images and darks in order of increasing exposure
;   (the resulting array lists the indexes of the exposures in increasing order).
imgsort = sort(imgsig[*,0])
dksort = sort(dksig[*,0])

for i = 0, (n/2)-1 do begin
  ; The indexes of the image and dark with the same exposure
  imgi = imgsort[i]
  dki = dksort[i]
  
  ; Subtract the dark (for each pixel) and add the resulting value to the signals array, which is now sorted
  for j = 0, (xsize*ysize)-1 do begin
    if (FLOAT(imgsig[imgi,j+1]) - FLOAT(dksig[dki,j+1])) lt 0 then begin
      signals[i,j] = 0
    endif else begin
      signals[i,j] = FLOAT(imgsig[imgi,j+1]) - FLOAT(dksig[dki,j+1])
    endelse
  endfor
  
  ; Add the corresponding exposure to the exposures array
  exposures[i] = imgsig[imgi,0]
endfor

;================================================================================

print, 'Plotting nonlinearity histograms...'
print, 'Please pick a directory to save distributions: '
newdir = dialog_pickfile(/directory)

READ, inc, PROMPT='Please provide the desired exposure interval between each nonlinearity histogram (ms). This should be a multiple of the exposure interval: '
; NOTE: Modify the interval value to whichever exposure interval the data was taken at.
interval = 0.05

; NOTE: Modify the starting and ending values of i depending on how many histograms you want to produce.
; For the NIRIS camera, the images reached saturation very early on in the testing process, so I shifted the range
;   for plotting histograms so that they would better represent the entire arc of the nonlinearity distributions.
; Typically, i would range from (n/4) to (n/2)-1 (roughly half saturation to full saturation).
for i = 10, (n/2)-10, FLOOR(inc/interval) do begin
  for j = 0, (xsize*ysize)-1 do begin
    ; Ignore pixels outside of the sun's area (pixels with nearly no increase in signal over the entire testing range).
    ; The threshold (1000 ADU) should be changed depending on the saturation value of the camera.
    if (avg(signals[*,j]) lt 1000) or (max(signals[*,j]) - min(signals[*,j]) lt 1000) then begin
      nonlins[j] = !VALUES.F_NAN
    endif else begin
      coeff = LINFIT(exposures[0:i], signals[0:i, j])
      mini = signals[i,j]
      maxi = 0.0
      maxsig = 0.0
      for index = 0, i do begin
        curr = signals[index,j] - (coeff[1]*exposures[index] + coeff[0])
        if curr gt maxi then begin
          maxi = curr
        endif
        if curr lt mini then begin
          mini = curr
        endif
        if signals[index,j] gt maxsig then begin
          maxsig = signals[index,j]
        endif
      endfor

      if maxsig eq 0 then begin
        nonlin = 100
      endif else begin
        nonlin = ((maxi+abs(mini))/maxsig)*100
      endelse
      nonlins[j] = nonlin
    endelse
  endfor

  ; Plot and save the distribution histogram to a specific directory.
  ; Ignore the pixels outside of the sun's area (where the nonlinearity is NaN).
  avgsig = 0.0
  count = 0UL
  for k = 0, xsize*ysize-1 do begin
    if finite(nonlins[k]) then begin
      avgsig += signals[i,k]
      count++
    endif
  endfor
  avgsig /= count

  mlin = 2*(median(nonlins))
  hist = HISTOGRAM(nonlins, LOCATIONS=binvals, BINSIZE=0.1, /NAN)
  histplot = PLOT(binvals, hist, XTITLE='Nonlinearity (%)', YTITLE='Frequency', $
    TITLE='Pixel Nonlinearity Distributions ('+strtrim(string(avgsig),2)+' ADU)', XRANGE=[0,mlin], /STAIRSTEP)
  savedir = newdir+'/nonlin_'+strtrim(string(avgsig),2)+'adu.jpg'
  histplot.Save, savedir
endfor

print, 'Finished.'

end