; Prompts the user to select a folder, and reads all FITS files into an array of file names
dir = dialog_pickfile(/directory,path='c:\')
list = file_search( dir + '*.fts' , count = n )

; Prompts the user to input the data's starting exposure and interval
READ, t0, PROMPT='Please input the starting exposure (in ms): '
READ, interval, PROMPT='Please input the exposure interval: '

exposures = ( FINDGEN(n/2) ) * interval + t0                
signals = fltarr(n/2) & darks = signals 

offset = 20

; If the dark field images are grouped at the end...
;for i = 0, (n/2)-1 do begin
;  img = readfits(list[i])
;  signals[i] = avg( img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 3] )
;endfor
;
;for j = n/2, n-1 do begin                                                                 
;  dk = readfits(list[j])                                                                       
;  darks[j] = avg( dk[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 3] )      
;endfor

print, 'Please input the coordinates of the center of the region of interest.'
READ, xcoord, PROMPT='x-coordinate: '
READ, ycoord, PROMPT='y-coordinate: '

; If the images and dark fields are taken together...
for i = 0, n/2 - 1 do begin
  img = readfits(list[2*i+1])
  dk = readfits(list[2*i])
  signals[i] = avg( img[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 3] )
  darks[i] = avg( dk[xcoord-offset : xcoord+offset, ycoord-offset : ycoord+offset, 3] )
endfor

res = signals - darks

; Plot the data and perform a linear regression
plot, exposures, res, psym=1, XTITLE='Exposure Time [ms]', YTITLE='Mean Signal [ADU]'

; Allows the user to test linear regression intervals until satisfied
cont = boolean(1)
while cont do begin
  READ, cil, PROMPT='Please input the linear max (in ms): '
  coeff = LINFIT(exposures[0:(cil*(1/interval))-1], res[0:(cil*(1/interval))-1])
  YFIT = coeff[0] + coeff[1]*exposures
  
  OPLOT, exposures, YFIT
  print, 'Equation of the linear regression from x =', t0, ' ms to x =', cil, ' ms: y =', coeff[1], 'x +', coeff[0]
  
  ; Determine the nonlinearity of the data
  mini = res[cil*(1/interval)-1]
  maxi = 0.0
  maxsig = 0.0
  for index = 0, (cil*(1/interval))-1 do begin
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
  
  print, 'The maximum deviation is ', maxi, ', the minimum deviation is ', mini, ' and the maximum signal is ', maxsig, '.'
  print, 'The nonlinearity over the linear range is ', nonlin*100, '%.'
  
  str = ''
  READ, str, PROMPT='Continue regression? (y/n) '
  if (str eq 'n') then begin
    cont = boolean(0)
  endif
endwhile

end

; possible improvement for the future: have the program read the first two FITS files and
; automatically set the initial exposure and interval