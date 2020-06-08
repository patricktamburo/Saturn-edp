;Data storage info
;
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/' 
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0277/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Titan/cors_0277/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'X'          ;Must be capitalized for comparison with data
fileID = 'xr94-63'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-ssa-rss-1-tboc5-v10/cors_0277/
if fileID eq 'kr173' then rsrfile = datapath+'S51TIOI2009173_1820NNNK55RD.1B2'
if fileID eq 'sr173' then rsrfile = datapath+'S51TIOI2009173_1820NNNS65RD.2B2'
if fileID eq 'xr173' then rsrfile = datapath+'S51TIOI2009173_1820NNNX55RD.1A2'
if fileID eq 'xl94-55' then rsrfile = datapath+'S49TIOC2009094_0233NNNX55LD.1B2'
if fileID eq 'xr94-55' then rsrfile = datapath+'S49TIOC2009094_0233NNNX55RD.1A2'
if fileID eq 'xl94-63' then rsrfile = datapath+'S49TIOC2009094_0233NNNX63LW.1N2'
if fileID eq 'xr94-63' then rsrfile = datapath+'S49TIOC2009094_0233NNNX63RW.1N2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/' 
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0278/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Titan/cors_0278/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'X'          ;Must be capitalized for comparison with data
fileID = 'xr173-14'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-ssa-rss-1-tboc5-v10/cors_0278/
if fileID eq 'xr173-65' then rsrfile = datapath+'S51TIOI2009173_1820NNNX65RD.2A2'
if fileID eq 'xr173-14' then rsrfile = datapath+'S51TIOI2009173_1845NNNX14RD.1A2'
if fileID eq 'xl173' then rsrfile = datapath+'S51TIOI2009173_1845NNNX14LD.2A2'
if fileID eq 'sr173' then rsrfile = datapath+'S51TIOI2009173_1845NNNS14RD.1B2'
if fileID eq 'sl173' then rsrfile = datapath+'S51TIOI2009173_1845NNNS14LD.2B2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/' 
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0276/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Titan/cors_0276/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'S'          ;Must be capitalized for comparison with data
fileID = 'sr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-ssa-rss-1-tboc5-v10/cors_0276/
if fileID eq 'xr' then rsrfile = datapath+'S49TIOC2009094_0233NNNX14RD.3A2'
if fileID eq 'sr' then rsrfile = datapath+'S49TIOC2009094_0233NNNS14RV.1N2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/' 
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0105/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Saturn/cors_0105/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'X'          ;Must be capitalized for comparison with data
fileID = 'xr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-s-rss-1-sroc-v10/cors_0105/
if fileID eq 'xr' then rsrfile = datapath+'s10sroe2005123_0740nnnx34rd.1a2'
if fileID eq 'sr' then rsrfile = datapath+'s10sroe2005123_0740nnns43rd.2b2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/'
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0106/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Saturn/cors_0106/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'S'          ;Must be capitalized for comparison with data
fileID = 'sr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-s-rss-1-sroc-v10/cors_0106/
if fileID eq 'sr' then rsrfile = datapath+'s10sroi2005123_0230nnns14rd.2b2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/'
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0107/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Saturn/cors_0107/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'X'          ;Must be capitalized for comparison with data
fileID = 'xr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-s-rss-1-sroc-v10/cors_0107/
if fileID eq 'xr' then rsrfile = datapath+'s10sroi2005123_0230nnnx14rd.2a2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/' 
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0203/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Saturn/cors_0203/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'X'          ;Must be capitalized for comparison with data
fileID = 'xr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-s-rss-1-sroc5-v10/cors_0203/
if fileID eq 'xr' then rsrfile = datapath+'s34saoe2007297_0745nnnx55rd.2a2'
if fileID eq 'sr' then rsrfile = datapath+'s34saoe2007297_0745nnns63rd.1b2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Titan/cors_0140/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'S'          ;Must be capitalized for comparison with data
fileID = 'sr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-ssa-rss-1-tboc1-v10/cors_0140/
if fileID eq 'sr' then rsrfile = datapath+'s19tioc2006078_0107nnns14rd.2a2'
if fileID eq 'xr' then rsrfile = datapath+'s19tioc2006078_0107nnnx14rd.1a2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;Paths for RSR files and figures/data output
if SCC_mode eq 0 then begin
  rootdatapath = '/Volumes/PW-2TB/'
  rootoutputpath = '/Users/paul/Documents/Boston_University/Research/'
endif
if SCC_mode eq 1 then begin
  rootdatapath = '/projectnb/cassini/'
  rootoutputpath = '/project/cassini/'
endif

datapath = rootdatapath+'Cassini_Ionosphere/Data/cors_0207/'
outputpath = rootoutputpath+'Cassini_Ionosphere/Output/Saturn/cors_0207/'
spawn, 'mkdir '+outputpath

;Choose band (x, s, or k) and circular polarization (l or r)
band = 'S'          ;Must be capitalized for comparison with data
fileID = 'sr'   ;A unique identifier for this file that matches an if condition below

;Choose the RSR file with the desired band and polarization
;The following data are from /pdsd/archive/data/co-ssa-rss-1-sroc5-v10/cors_0207/
if fileID eq 'sr' then rsrfile = datapath+'s36saoe2007353_0524nnns63rd.1b2'
if fileID eq 'xr-55' then rsrfile = datapath+'s36saoe2007353_0524nnnx55rd.2a2'
if fileID eq 'xr-63' then rsrfile = datapath+'s36saoe2007353_0524nnnx63rd.1a2'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

