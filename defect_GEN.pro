darkname = dialog_pickfile()
dark = readfits(darkname)

dev = stddev(dark)
sig = avg(dark)

window, 0, xsize=512, ysize=512
dark0 = rebin(dark, 512, 512)
tvscl, dark(*,*,0)

;print, dev, sig
warm = where(dark GT sig + 5.*dev)
hot = where(dark GT 10000.)

print, 'warm pixels: ', warm
print, 'hot pixels: ', hot

end