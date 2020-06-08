pro vline, x, y, label=label, pos=pos, _extra = _extra

; plots a vertical line(s) at the position(s) specified, on the current graph
; works transparentally with logarithmic plots.
;
; INPUT Variable
;	x	: the location(s) on the horizontal (x) axis to place a line [scalar or vector]
; 	y 	: [optional] a 2-element array specifying the start and stop y locations of the lines,
;			in data coordinates.
; OPTIONAL KEYWORDS
;	label 	:	String label(s) for the line, to appear to the LEFT of the line [scalar or vector].
;	pos 	:	The normalized vertical location of the label(s) [default = 0.5].
;	_extra	:	Any graphics keywords accepted by either plots OR xyouts (COLOR, THICK, CHARTHICK, etc)
;
; RETURN VALUE: none
;

nx = n_elements(x)
if n_elements(y) eq 0 then $
	if !y.type eq 1  then y = [10^(!y.crange)] else y = !y.crange
if !x.type then xcr = [10^(!x.crange)] else xcr = !x.crange

for i = 0, nx-1 do begin
	if (x[i] LT xcr[1]) AND (x[i] GT xcr[0]) then begin
		plots, x[i], y[0], _extra= _extra, noclip=0
		plots, x[i], y[1], /continue, _extra = _extra, noclip=0
	endif
endfor

if keyword_set(label) then begin
	if n_elements(pos) eq 0 then pos = fltarr(nx)+0.5
	if n_elements(pos) LT nx then pos_ = replicate(pos[0], nx) else pos_=pos
	if n_elements(label) LT nx then label_ = replicate(label[0],nx) else label_=label
	for i=0, nx-1 do begin
		; pos should contain the normalized coord of the y-position of the label
		; will appear in same color as line (to the left of the line)
		if !x.type then f1 = alog10(x[i]) else f1 = x[i]
		lp = (f1 - xcr[0])/(xcr[1] - xcr[0])
		lp = lp - 0.01 ; a little to the left
		lp = xcr[0] + lp * (xcr[1] - xcr[0])
		if !x.type then lpx = 10^lp else lpx = lp
		lpy = !y.crange[0] + pos_[i]*(!y.crange[1]-!y.crange[0])
		if !y.type then lpy = 10^lpy
		xyouts, lpx, lpy, align=1.0, label_[i], _extra=_extra
	endfor
endif

end

