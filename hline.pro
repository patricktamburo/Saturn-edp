pro hline, y, x, label=label, pos=pos, _extra = _extra

; plots a horizontal line(s) at the position(s) specified, on the current graph
; works transparentally with logarithmic plots.
; can also label the line or lines.
;
; INPUT Variable
;	y	: the location(s) on the vertical (y) axis to place a line [scalar or vector]
;   x   : [optional] a 2-element array specifying the start and stop x locations of the lines,
;			in data coordinates.
; OPTIONAL KEYWORDS
;	label 	:	String label(s) for the line, to appear ABOVE the line [scalar or vector].
;	pos 	:	The normalized horizontal location of the label(s) [default = 0.05].
;	_extra	:	Any graphics keywords accepted by either plots OR xyouts (COLOR, THICK, CHARTHICK, etc)
;
; RETURN VALUE: none
;
; EXAMPLE:	Say you have a graph and want to plot 2 horizontal lines on it, at values y=2 and y=6.
;			You also want to label them with labels "Line 1" and "Line 2", and make the lines be thick,
;			dashed, and red.
;
;		hline, [2,6], label=['Line 1','Line 2'], linestyle=2, color = 150, thick=3
;
;  If you are motivated, make it so the keywords to _extra can be vectors as well, modify the program to do
;  this and let me know!
;
;		Chris
;
ny = n_elements(y)
if n_elements(x) eq 0 then $
	if !x.type then x = [10^(!x.crange)] else x = !x.crange

if !y.type then ycr = [10^(!y.crange)] else ycr = !y.crange

for i = 0, ny-1 do begin
	if (y[i] LT ycr[1]) AND (y[i] GT ycr[0]) then begin
		plots, x[0], y[i], _extra= _extra, noclip=0
		plots, x[1], y[i], /continue, _extra = _extra,noclip=0
	endif
endfor

if keyword_set(label) then begin
	if n_elements(pos) eq 0 then pos = fltarr(ny)+0.05
	if n_elements(pos) LT ny then pos_ = replicate(pos[0], ny) else pos_=pos
	if n_elements(label) LT ny then label_ = replicate(label[0],ny) else label_=label
	for i=0, ny-1 do begin
		; pos should contain the normalized coord of the y-position of the label
		; will appear in same color as line (to the left of the line)
		if !y.type then f1 = alog10(y[i]) else f1 = y[i]
		lp = (f1 - ycr[0])/(ycr[1] - ycr[0])
		lp = lp + 0.02 ; a little above
		lp = ycr[0] + lp * (ycr[1] - ycr[0])
		if !y.type then lpy = 10^lp else lpy = lp
		lpx = !x.crange[0] + pos_[i]*(!x.crange[1]-!x.crange[0])
		if !x.type then lpx = 10^lpx
		xyouts, lpx, lpy, align=0.0, label_[i], _extra=_extra
	endfor
endif

end

