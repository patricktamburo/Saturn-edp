pro rsr_proc_single

;+
;Name:
;     rsr_proc_single
;Purpose:
;     This program reads RSR data files from Cassini RSS
;     and extracts frequency vs. time for the observations.
;     It is meant to be run on specific jobs, not batches.
;Calling Sequence:
;     rsr_proc_single
;Inputs:
;     Location of the data files, band, polarization, and 
;     other choices about the analysis method.
;Output:
;     Various figures showing the I/Q values, power spectra,
;     and the extracted frequency vs. time relation. Also an
;     IDL save file containing the frequency vs. time result.
;Keywords:
;     None
;Authors:
;     Paul Withers, Paul A. Dalba
;     Astronomy Department, Boston University
;     January 2017
;-

; Record time program starting running
spawn, 'date', codestarttime

;----------------------------------------START USER INPUTS----------------------------------------;
;Is this going to be run on the SCC? 1=True, 0=False
SCC_mode = 0

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/'
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0140/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Saturn/cors_0140/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'X'          ;Must be capitalized for comparison with data
fileID = 'xr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-ssa-rss-1-sroc5-v10/cors_0207/
if fileID eq 'xr' then rsrfile = datapath+'s19tioc2006078_0107nnnx14rd.1a2'

;Set the number of I/Q 'sample words' in the RSR data file. This can be found in the corresponding LBL file.
; For now, hold this at 4000L and stop the code if nwords does NOT equal 4000L
nwords = 4000L

;Set the number of samples that will be operated upon in each FFT. This divided by the sampling rate determines
; the duration of data in each FFT.
nforfft = 16000L

;Set the number of zeros required for zero padding. Computational efficient padding uses a power of 2.
npad = 2LL^25LL

;Determine the rows (SFDUs) of the RSR file to analyze in total. This is equivalent to setting the time bounds
; of the returned frequency vs. time.
; The (iend-istart) should be an integer multiple of (nforfft/nwords). This quotient is the number of rows that
; go into every FFT, so the total number of rows should be an integer multiple of this.
; Enter istart = 'AUTO' if you would like the code to automatically determine the istart/iend values to sample the entire dataset. 
; Otherwise, enter your own values as described above.
istart = 'auto'
;istart = 'auto'   
iend = 52000  
;----------------------------------------END USER INPUTS----------------------------------------;

;Start just by printing the RSR file this run is for. This is useful for SCC record keeping
print, 'This run is for RSR file:'
print, rsrfile

;Open the corresponding label file to extract the number of rows and check the number of words.
; Extract the file prefix assuming typical RSR naming convention
fileprefix = (strsplit((strsplit(rsrfile,'/',/extract))(-1),'.',/extract))(0)
;Open the label file. Must be in the same directory as the RSR file. (Annoyingly) text case can vary.
; This means the label file can be .LBL or .lbl. Open the lable by matching the case of the prefix.
if fileprefix eq strupcase(fileprefix) then $
  openr, lun0, datapath+fileprefix+'.LBL', /get_lun $
else $
  openr, lun0, datapath+fileprefix+'.lbl', /get_lun
line = ''
while not EOF(lun0) do begin
  readf, lun0, line
  ;Search for the number of rows
  if strmatch(line,'*ROWS*') then nrows = long((strsplit(line,'=',/extract))(-1))
  ;Search for the number of sample words
  if strmatch(line,'*SAMPLE WORDS*') then begin
    readf, lun0, line
    readf, lun0, line
    readf, lun0, line
    samplewords = long((strsplit(line,'=',/extract))(-1))
  endif
endwhile
free_lun, lun0

stop

;Run a check to make sure the nwords=4000 as expected
if samplewords ne nwords then begin
  print, 'ERR6: Number of sample words not equal to 4000.'
  stop
endif

;Run a check to make sure this is a one-way occ. Stop if two-way.
if not strmatch(fileprefix, '*NNN*') then begin
  if not strmatch(fileprefix, '*nnn*') then begin
    print, 'ERR7: This RSR file appears to be for a two-way occulation.'
    stop
  endif  
endif

;With nrows, automatically determine the istart/iend values if the user so desires
if string(istart) eq 'auto' then begin
  istart = 0L
  iend = floor(nrows*nwords/nforfft)*nforfft/nwords
endif


;Define arrays for storage of data read from binary RSR file
majorclassarr = intarr(nrows)  ;Major data class
minorclassarr = intarr(nrows)  ;Minor data class

rsnarr = intarr(nrows)          ;Record sequence number 
sctarr = intarr(nrows)          ;Spacecraft ID
proccenterarr = intarr(nrows)   ;Signal processing center ID
rsrreceiverarr = intarr(nrows)  ;Radio receiver ID
upfreqarr = strarr(nrows)       ;Uplink frequency band
downfreqarr = strarr(nrows)     ;Downlink frequency band

sfdusampleratearr = intarr(nrows)  ;Sample rate (kHz)
rftoifmhzarr = intarr(nrows)       ;RF to IF down conversion frequency (MHz)
ddclomhzarr = intarr(nrows)        ;DDC LO frequency (MHz)
dataerrorarr = intarr(nrows)       ;Data error count

digadcyear = intarr(nrows)  ;UTC year of ADC
digadcdoy = intarr(nrows)   ;UTC day of year of ADC
digadcsec = dblarr(nrows)   ;UTC second of data of ADC
sfduyear = intarr(nrows)    ;UTC year of SFDU data and models
sfdudoy = intarr(nrows)     ;UTC day of year of SFDU data and models
sfdusec = dblarr(nrows)     ;UTC second of day of SFDU data and models

f1arr = dblarr(nrows)  ;Sub-channel frequency polynomial coeff. 1 (Hz)
f2arr = dblarr(nrows)  ;Sub-channel frequency polynomial coeff. 2 (Hz)
f3arr = dblarr(nrows)  ;Sub-channel frequency polynomial coeff. 3 (Hz)

scfreq1 = dblarr(nrows)  ;Sub-channel frequency (Hz) at beginning of second
scfreq2 = dblarr(nrows)  ;Sub-channel frequency (Hz) at middle of second
scfreq3 = dblarr(nrows)  ;Sub-channel frequency (Hz) at end of second

rsrqarr = intarr(nrows, nwords)   ;All Q values from all rows
rsriarr = intarr(nrows, nwords)   ;All I values from all rows
rsrindex = indgen(nwords)*2L      ;Even indices of rsr Q/I values


;Open the data file to read
openr, lun1, rsrfile, /get_lun

;Begin a loop that will define each column as a variable and read in the column data
row=0L
while row lt nrows do begin
  col01 = '    '
  col02 = ' '
  col03 = ' '
  col04 = 0S
  col05 = '    '
  col06 = 0L
  col07 = 0L
  col08 = 0S
  col09 = 0S
  col10 = 0S
  col11 = 0S
  col12 = 0B  ;Major data class. Should always be set to '21' for radio science data
  col13 = 0B  ;Minor data class. Should always be set to '4' for radio science RSR data
  col14 = 0B
  col15 = 0B
  col16 = 0S
  col17 = 0S
  col18 = 0B
  col19 = 0B
  col20 = 0S
  col21 = 0S  ;Record sequence number. Should increment by 1 for each SFDU (i.e. row)
  col22 = 0B  ;Signal processing center (i.e. DSN site)
  col23 = 0B
  col24 = 0B  ;Radio science receiver (i.e. which receiver and rack were used). Combine with col22 to uniquely identify hardware.
  col25 = 0B
  col26 = 0B
  col27 = 0B  ;Spacecraft identifier
  col28 = 0S
  col29 = ' ' ;Uplink frequency. Could be 'S', 'X', or 'K'
  col30 = ' ' ;Downlink frequency. Could be 'S', 'X', or 'K'
  col31 = 0B
  col32 = 0B
  col33 = 0B
  col34 = 0B
  col35 = 0B
  col36 = 0B
  col37 = 0B
  col38 = 0B
  col39 = 0S  ;UTC year of ADC 
  col40 = 0S  ;UTC day of year of ADC 
  col41 = 0L  ;UTC second of day of ADC 
  col42 = 0B
  col43 = 0B  ;Data error count. Anything >0 may suggest data are corrupted by hardware error. 
  col44 = 0S  ;Sample rate in kilosample per second
  col45 = 0S  ;DDC LO frequency (MHz). Needed to reconstruct sky frequency
  col46 = 0S  ;RF to IF DDC LO frequency (MHz). Needed to reconstruct sky frequency
  col47 = 0S  ;UTC year of SFDU data and models
  col48 = 0S  ;UTC day of year of SFDU data and models
  col49 = 0D  ;UTC second of day of SFDU data and models
  col50 = 0D
  col51 = 0D
  col52 = 0D
  col53 = 0D
  col54 = 0D
  col55 = 0D  ;Radio frequency (Hz) at beginning of second
  col56 = 0D  ;Radio frequency (Hz) at middle of second
  col57 = 0D  ;Radio frequency (Hz) at end of second
  col58 = 0D  ;Sub-channel frequency (Hz) at beginning of second
  col59 = 0D  ;Sub-channel frequency (Hz) at middle of second
  col60 = 0D  ;Sub-channel frequency (Hz) at end of second
  col61 = 0D  ;Sub-channel frequency polynomial coeff. 1 (Hz)
  col62 = 0D  ;Sub-channel frequency polynomial coeff. 2 (Hz)
  col63 = 0D  ;Sub-channel frequency polynomial coeff. 3 (Hz)
  col64 = 0D
  col65 = 0D
  col66 = 0D
  col67 = 0D
  col68 = 0D
  col69 = bytarr(16)
  col70 = 0S
  col71 = 0S
  col72 = intarr(nwords*2)  ;Quadrature (Q) and in-phase (I) data. Each has nwords of data
  ;Now read the binary data in order
  readu, lun1, col01
  readu, lun1, col02
  readu, lun1, col03
  readu, lun1, col04
  readu, lun1, col05
  readu, lun1, col06
  readu, lun1, col07
  readu, lun1, col08
  readu, lun1, col09
  readu, lun1, col10
  readu, lun1, col11
  readu, lun1, col12
  readu, lun1, col13
  readu, lun1, col14
  readu, lun1, col15
  readu, lun1, col16
  readu, lun1, col17
  readu, lun1, col18
  readu, lun1, col19
  readu, lun1, col20
  readu, lun1, col21
  readu, lun1, col22
  readu, lun1, col23
  readu, lun1, col24
  readu, lun1, col25
  readu, lun1, col26
  readu, lun1, col27
  readu, lun1, col28
  readu, lun1, col29
  readu, lun1, col30
  readu, lun1, col31
  readu, lun1, col32
  readu, lun1, col33
  readu, lun1, col34
  readu, lun1, col35
  readu, lun1, col36
  readu, lun1, col37
  readu, lun1, col38
  readu, lun1, col39
  readu, lun1, col40
  readu, lun1, col41
  readu, lun1, col42
  readu, lun1, col43
  readu, lun1, col44
  readu, lun1, col45
  readu, lun1, col46
  readu, lun1, col47
  readu, lun1, col48
  readu, lun1, col49
  readu, lun1, col50
  readu, lun1, col51
  readu, lun1, col52
  readu, lun1, col53
  readu, lun1, col54
  readu, lun1, col55
  readu, lun1, col56
  readu, lun1, col57
  readu, lun1, col58
  readu, lun1, col59
  readu, lun1, col60
  readu, lun1, col61
  readu, lun1, col62
  readu, lun1, col63
  readu, lun1, col64
  readu, lun1, col65
  readu, lun1, col66
  readu, lun1, col67
  readu, lun1, col68
  readu, lun1, col69
  readu, lun1, col70
  readu, lun1, col71
  readu, lun1, col72
  
  ;Store column data in arrays and swap endian
  majorclassarr(row) = swap_endian(col12)  ;Major data class. Should always be set to '21' for radio science data
  minorclassarr(row) = swap_endian(col13)  ;Minor data class. Should always be set to '4' for radio science RSR data
  
  rsnarr(row) = swap_endian(col21)          ;Record sequence number
  sctarr(row) = swap_endian(col27)          ;Spacecraft identifier
  proccenterarr(row) = swap_endian(col22)   ;Signal processing center ID
  rsrreceiverarr(row) = swap_endian(col24)  ;Radio receiver ID
  upfreqarr(row) = swap_endian(col29)       ;Uplink frequency band
  downfreqarr(row) = swap_endian(col30)     ;Downlink frequency band
  
  sfdusampleratearr(row) = swap_endian(col44)  ;Sample rate (kilosamples/sec) of this SFDU (i.e. row)
  rftoifmhzarr(row) = swap_endian(col46)       ;RF to IF down conversion frequency in MHz
  ddclomhzarr(row) = swap_endian(col45)        ;DDC LO frequency in MHz
  dataerrorarr(row) = swap_endian(col43)       ;Data error count
  
  digadcyear(row) = swap_endian(col39)  ;UTC year of ADC
  digadcdoy(row) = swap_endian(col40)   ;UTC day of year of ADC
  digadcsec(row) = swap_endian(col41)   ;UTC second of day of ADC
  sfduyear(row) = swap_endian(col47)    ;UTC year of SFDU data and models
  sfdudoy(row) = swap_endian(col48)     ;UTC day of year of SFDU data and models
  sfdusec(row) = swap_endian(col49)     ;UTC second of day of SFDU data and models
  
  f1arr(row) = swap_endian(col61)   ;Sub-channel frequency polynomial coeff. 1 (Hz)
  f2arr(row) = swap_endian(col62)   ;Sub-channel frequency polynomial coeff. 2 (Hz)
  f3arr(row) = swap_endian(col63)   ;Sub-channel frequency polynomial coeff. 3 (Hz)
  
  scfreq1(row) = swap_endian(col58)   ;Sub-channel frequency (Hz) at beginning of second
  scfreq2(row) = swap_endian(col59)   ;Sub-channel frequency (Hz) at middle of second
  scfreq3(row) = swap_endian(col60)   ;Sub-channel frequency (Hz) at end of second
  
  rsrq = swap_endian(col72(rsrindex))  ;Select Q values, which are in the even indices
  rsri = swap_endian(col72(rsrindex+1L))  ;Select I values, which are in the odd indices
  
  ;Add the Q and I data from this row to their master arrays 
  rsrqarr(row,*) = rsrq(*)
  rsriarr(row,*) = rsri(*)
  
  row++
endwhile

;Release this unit
free_lun, lun1

;Save the I/Q pairs and the sfdusec for this particular RSR
;save, filename='/Users/paul/Desktop/IQ_pairs.sav', rsrqarr, rsriarr, sfdusec
;stop

;Complete multiple checks of the data and user inputs before continuing analysis
;First, check to make sure the sample rate for this block of data did not change.  
if n_elements(sfdusampleratearr) ne n_elements(where(sfdusampleratearr eq sfdusampleratearr[0])) then begin 
  print, 'ERR1: Sample rate not constant.'
  stop
 endif 

;Second, check that the band is as expected, and that it does not change
if n_elements(sfdusampleratearr) ne n_elements(where(downfreqarr eq band)) then begin
  print, 'ERR2: Downlink frequency mismatch.'
  stop
endif

;Third, check that the data error count is zero. Anything more could indicate corrupted data
if max(dataerrorarr) ne 0 then begin
  print, 'ERR3: Data error count is nonzero.'
  stop
endif

;Fourth, check that the user-defined values will result in even FFT blocks
if (nforfft mod nwords) ne 0 then begin
  print, 'ERR4: nforfft not chosen correctly.'
  stop
endif
if ((iend-istart) mod (nforfft/nwords)) ne 0 then begin
  print, 'ERR5: i-range not chosen correctly.'
  stop
endif


;Define the sample rate from the array, assumign the sample rate check has been passed. Note the unit of kHz.
samplerate = sfdusampleratearr[0]*1000D

;t is an npad-sized array from 0 to (npad/sample rate). It's the times matching the data to be operated on.
t = dindgen(npad)/(npad-1LL) * npad / samplerate
trange = max(t) - min(t)
dt = t(1) - t(0)

;Just a subset of the full time range - matching nforfft - and call this subset torig. This is the time 
; range at the proper sampling for just the data without the padded zeros.
torig = t[0:nforfft-1]
torigrange = max(torig) - min(torig)
dtorig = torig[1] - torig[0]

;Create the frequency arrays for both cases (full t and subset t). k counts up from zero to a positive max, 
; then reverses down to zero with negative numbers.
; Given the name k, although not a wavenumber, to distinguish between all the frequencies in this code.
k = [l64indgen(npad/2LL + 1LL), -1LL*reverse(1LL + l64indgen(npad/2LL - 1LL))]/trange
kdummy = [l64indgen(nforfft/2LL + 1LL), -1LL*reverse(1LL + l64indgen(nforfft/2LL - 1LL))]/torigrange


;Define arrays that will be the output products of the FFTs. The size is the number of FFTs to be run: nout
; delta-i is the number of sfdu's, times by nwords per sfdu gives total words, then dividing by the number 
; of samples for the FFT gives the number of FFTs we have to do to get through all i rows of data
nout = 1D * (iend-istart) * nwords / nforfft  

;Timing information
sfduyearout = intarr(nout)
sfdudoyout = intarr(nout)
sfdusecout = dblarr(nout)

;Information needed to upsample the frequency later
rftoifmhzout = intarr(nout)
ddclomhzout = intarr(nout)
ncofreqout = dblarr(nout)

;The output frequency from the FFT
kposaaout = dblarr(nout)
;The maximum power in the power spectrum. If small, then the observations only sampled noise. If large, the observations sampled something physical.
igcomplexaaout = dcomplexarr(nout)


;Begin a loop over kindex of size nout to go through the selected data and perfrom FFTs
i=istart 
kindex=0L ; counter within this loop

while i lt iend do begin
  print, i, kindex, nout
   
  ;Since nout is the number of FFT operations that must be done, the nblock is the number of rows (SFDUs) in each of these nout runs
  ; This is equal to nforfft/nowrds, which we know is by definition an integer
  nblock = (iend - istart) / nout
  
  ;Since the kindex will only be evaluated every nblock i-steps, fill these ancillary data arrays with the midpoint values 
  sfduyearout[kindex] = sfduyear(i+nblock/2L)
  sfdudoyout[kindex] = sfdudoy(i+nblock/2L)
  sfdusecout[kindex] = sfdusec(i+nblock/2L)
  rftoifmhzout[kindex] = rftoifmhzarr(i+nblock/2L)
  ddclomhzout[kindex] = ddclomhzarr(i+nblock/2L)
  ;The subchannel frequency polynomial evaluates frequency over 1 second intervals. tms is the number of milliseconds ellapsed since the last whole second.
  tms = 1D3*(sfdusecout[kindex]-floor(sfdusecout[kindex]))
  thisf1 = f1arr(i+nblock/2L)
  thisf2 = f2arr(i+nblock/2L) ; calculate F1, F2, F3 at X.000 seconds, so F2 and F3 irrelevant
  thisf3 = f3arr(i+nblock/2L)
  ;Calculate the final frequency for this block
  ncofreqout[kindex] = thisf1 + thisf2*((tms+0.5)/1D3) + thisf3*((tms+0.5)/1D3)^2
  
  ;Define I and Q value arrays
  ivorig = dblarr(nforfft)
  qvorig = dblarr(nforfft)
  
  ;The following loop selects the I/Q values from the original full array, and reshapes them into ivorig and qvorig. 
  block=0
  while block lt nblock do begin
    ivorig[block*nwords:(block+1L)*nwords - 1L] = rsriarr(i+block,*) 
    qvorig[block*nwords:(block+1L)*nwords - 1L] = rsrqarr(i+block,*)
    block++
  endwhile
  
  ;Set up iv and ivcomplex with the zero padding
  iv = dblarr(npad)
  ivcomplex = dcomplexarr(npad)
  ;Fill the ivcomplex array with the I/Q values, which are the real and complex measurements
  ivcomplexorig = complex(ivorig, qvorig)
  ;Fill the first nforfft elements of ivcomplex with the original I/Q values. The rest will be left as padding zeros.
  ivcomplex[0:nforfft-1] = ivcomplexorig
  
  ;Perform the FFT on the complex array with the nforfft I/Q values and npad worth of padded zeros
  igcomplex = fft(ivcomplex)
  ;Run another FFT but only on the nforfft data, without the zero padding.
  igcomplexorig = fft(ivcomplexorig)
  
  ;First the element where the power spectrum was maximum (aa?)
  aa = where(abs(igcomplex) eq max(abs(igcomplex)))
  aa = aa[0]
  
  ;With the index of the max power from above, determine the maximum frequency k. 
  kposaaout[kindex] = k[aa]
  ;With the index of the max power, also save just what that max power is. This max power will still be low if the observations only sampled background signal.
  igcomplexaaout[kindex] = igcomplex[aa]
  
  
  ;Plot the original nforfft I/Q values. Turn this off if running all FFTs or IDL will break.
  ;window, 0
  ;plot, torig, ivorig
  ;oplot, torig, qvorig, color=255
  
  ;Additional plots of a portion of the I/Q values and the two power spectra (w/ and w/o zero padding).
  plotit1=0
  if SCC_mode eq 1 then plotit1 = 0
  if plotit1 eq 1 then begin
    set_plot, 'ps'
     
    !p.position=0
    !p.charsize=1.5
    !p.charthick=2.0
    !p.font=0
    !p.thick=5.0
    !x.thick=5.0
    !y.thick=5.0
    !z.thick=5.0
     
    tvlct, [0,255,0,0,150,255,0,255,255,128,128], [0,0,255,0,150,255,255,0,255,128,0],$
     [0,0,0,255,150,255,255,255,0,0,128]
    
    device, /portrait, encapsulated=0, color=1, filename=strcompress(outputpath+'iqsample_'+strtrim(string(kindex),1)+'.ps'), /bold
    plot, torig, ivorig, xra=[0,0.5], psym=1, symsize=0.2, xtitle='Time (seconds)', ytitle='I and Q'
    oplot, torig, qvorig, color=4, psym=1, symsize=0.2
    xyouts, 0.05, 1.5e4, 'I=Black, Q=Grey'
    device, /close
    
    device, /portrait, encapsulated=0, color=1, filename=strcompress(outputpath+'fftsample_'+strtrim(string(kindex),1)+'.ps'), /bold
    plot, k, abs(igcomplex)^2/max(abs(igcomplex)^2), /ylog, xra=[6,12], yra=[1e-8,1e0], xtitle='Frequency (Hz)', ytitle='Normalized Power'
    ;Plot the nforftt only FFT using kdummy, which has nforfft points
    oplot, kdummy, abs(igcomplexorig)^2/max(abs(igcomplexorig)^2), color=4, psym=-1
    xyouts, 7.25, 1e-7, 'Grey = Without zero padding', color=4
    xyouts, 7.25, 3e-8, 'Black = With zero padding'
    device, /close
    
    if !d.name eq 'PS' then device, /close
    !p.position=0
    !p.multi=0
    set_plot, 'x'
    !p.thick=1
    !x.thick=1
    !y.thick=1
    !x.margin=[10,3]
    !y.margin=[4,2]
    !p.color=-1  
  endif
  
  ;Increase the kindex to step to the next block of data and determine a new time and frequency
  kindex++
  ;Increase the i index. It should jump forward the number of rows in one block.
  i = istart + nblock * kindex
endwhile


;Set up a directory in which to save the output of this run
rsrdir = strsplit(rsrfile, '/', /extract) 
spawn, 'mkdir '+outputpath+rsrdir[-1]
savepath = outputpath+rsrdir[-1]+'/'

;Change the yrange for the freq. plot based on the band
if band eq 'X' then ylim=[-40,40]
if band eq 'S' then ylim= [-15,15]
if band eq 'K' then ylim=[20,50]

;Plot the final frequency vs. time of the portion of data analyzed
if SCC_mode eq 0 then begin
  window, 1
  plot, sfdusecout, kposaaout, yra=ylim
  
  ;Plot the max power for each FFT as a function of time
  window, 2
  plot, sfdusecout, abs(igcomplexaaout)^2
endif

;Save these plots
plotit2 = 1
if SCC_mode eq 1 then plotit2 = 0
if plotit2 eq 1 then begin
  set_plot, 'ps'
  
  !p.position=0
  !p.charsize=1.5
  !p.charthick=2.0
  !p.font=0
  !p.thick=5.0
  !x.thick=5.0
  !y.thick=5.0
  !z.thick=5.0
  
  device, /landscape, encapsulated=0, color=1, filename=savepath+'freq_vs_time.ps', /bold
  plot, sfdusecout, kposaaout, xtitle='SFDU Time (seconds)', ytitle='Frequency (Hz)'
  device, /close
  
  device, /landscape, encapsulated=0, color=1, filename=savepath+'max_power_vs_time.ps', /bold
  plot, sfdusecout, abs(igcomplexaaout)^2, /ylog, xtitle='SFDU Time (seconds)', ytitle='Max Power'
  device, /close
  
  !p.position=0
  !p.multi=0
  set_plot, 'x'
  !p.thick=1
  !x.thick=1
  !y.thick=1
  !x.margin=[10,3]
  !y.margin=[4,2]
  !p.color=-1
endif

;Also the important variables in a save file
save, filename=savepath+'output.sav', sfduyearout, sfdudoyout, sfdusecout, rftoifmhzout, ddclomhzout, ncofreqout, kposaaout, igcomplexaaout

;If this is a full run (1-sec sampling with 2^25 zeros), then make a second directory and copy the output file. 
; This protects from accidental overwrite
if nforfft/samplerate eq 1 then begin
  if npad eq 2LL^25LL then begin
    print, 'Duplicating full output in separate directory as well'
    spawn, 'mkdir '+savepath+'full_output'
    save, filename=savepath+'full_output/output.sav', sfduyearout, sfdudoyout, sfdusecout, rftoifmhzout, ddclomhzout, ncofreqout, kposaaout, igcomplexaaout
  endif  
endif

spawn, 'date', codeendtime
print, ''
print, codestarttime
print, codeendtime

if SCC_mode eq 0 then begin
  stop
endif
 
end
