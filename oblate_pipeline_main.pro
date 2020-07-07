pro oblate_pipeline_main
  ;Author
  ;   Pat Tamburo, Boston University, June 2020
  ;Purpose
  ;   Runs simple_occ and freqwork on processed RSR data, utilizing an oblate representation of Saturn. 
  
  body =     'Saturn'    ;Options are 'Titan' or 'Saturn' at present
  cors = ['106']
  flyby =    'S07'
  obs =      'S07N'
  bands_array =    ['SX', 'SX', 'XK', 'XK', 'XK']
  station_array =  ['14', '43', '25', '26', '34']
  
  print, 'Beginning simple_occ itloop = 0 for cors ',cors
  oblate_pipeline_simple_occ, 0, cors, body
  print, ''
  
  for i = 0, n_elements(bands_array)-1 do begin
    bands = bands_array[i]
    station = station_array[i]
    print, 'Beginning freqwork itloop = a for ',obs,' ', bands,' ',station
    oblate_pipeline_freqwork, 'a', body, flyby, obs, bands, station
    print, ''
  endfor
  
  print, 'Beginning simple_occ itloop = 1 for cors ',cors
  oblate_pipeline_simple_occ, 1, cors, body
  print, ''
  
  print, 'Beginning simple_occ itloop = 2 for cors ',cors
  oblate_pipeline_simple_occ, 2, cors, body
  print, ''

  for i = 0, n_elements(bands_array)-1 do begin
    bands = bands_array[i]
    station = station_array[i]
    print, 'Beginning freqwork itloop = b for ',obs,' ', bands,' ',station
    oblate_pipeline_freqwork, 'b', body, flyby, obs, bands, station
    print, ''
  endfor
  
  
stop
end