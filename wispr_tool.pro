;;
;
; IDL tool to generate WISPR orbital ephemeris with relevant parameters
; for tomography. It also generates synthetic headers for the
; corresponding fits files. by A.M.Vasquez (albert@iafe.uba.ar)
;
; INPUTS:  work in progress.
;
; OUTPUTS: work in progress.
;  
; CALL SEQUENCE: wispr_tool,[/option1,/option2,....]
;
; HISTORY
; A.M. Vasquez   - Oct-15-2017 - Version 1.0
;                - Oct-31-2017 - Version 2.0, first fully functional version,
;                  that produces Blank FITS files for selected dates.
;                - Nov-17-2017 - Version 3.0. Accurate computation of
;                  sun-fov-size in P.O.S. (v2.0 has an approximation),
;                  also corrected computation of CRPIX1.
;                - Nov-19-2017 - Version 3.1. Added /squareFOV so that
;                  square FOVs are now optional. Adapted needed code
;                  to properly function with both square and non-square cases.
;                - Dec-11-2017 - Version 3.2. Added /bin option to
;                  generate images of smaller size, rebined by "bf",
;                  default value is bf=4.
;                - Feb-05-2018 - Version 4.0. Adding non-constant
;                  cadence to achieve uniform longitudinal coverage.
;
;;

pro wispr_tool,loadk=loadk,correction=correction,FullList=FullList,ShortList=ShortList,$
  SciOrbBrief=SciOrbBrief,$
  ExtendedOrbits=ExtendedOrbits,ScienceOrbits=ScienceOrbits,ConstantCadence=ConstantCadence,Cadence=Cadence,UniformLong=UniformLong,DeltaLong=DeltaLong,$
  CircularOrbits=CircularOrbits,Equatorial=Equatorial,OffEquator=OffEquator,$
  CreateFITS=CreateFITS,Outdir=Outdir,basedir=basedir,$
  SquareFOV=SquareFOV,bin=bin,bf=bf

common constants,c,rsun,au
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,px_inner_arcsec,px_outer_arcsec,sun_spp_vector_HCI,sun_spp_vector_HAE,sun_spp_vector_HEE,sun_spp_vector_HEQ
common output,listtype
common indexes,i,j
common SynthFITS,hdr_Inner_0,hdr_Outer_0,img_Inner_0,img_Outer_0,datadir,epoch,et
common FOVs,FOVW_inner_px,FOVH_inner_px,FOVW_outer_px,FOVH_outer_px
common suffixes,suffix_binfac,suffix_square
common experiment,CircularOrbit_flag,CircularOrbit_Radius,CircularOrbit_Carlon,CircularOrbit_Carlat,experiment_suffix
common FOV_POINTINGS,alphaI_C,deltaI,alphaI_E,alphaI_W,alphaO_C,deltaO,alphaO_E,alphaO_W ; all in deg

  if NOT keyword_set(Outdir) then Outdir='../TestImage/'
  ;Handle basedir and datadir
  if not keyword_set(basedir) then basedir='/data1/'
  datadir=Outdir

  ;load WISPR FOVs
   load_fov_pointings
  
  ;load useful constants
   loadconstants

  ;LOAD KERNELS
  if keyword_set(loadk) then begin
      ;; Load kernels to initialize program.
      ;; Update by W. Thompson March 2017 to leverage new sunspice.
      load_sunspice_gen  
      ;; Load a set of kernels:  Use a meta kernel for convenience.
      cspice_furnsh, 'wispr_albert.tm'
   endif

 ;Load sample FITS
  mreadfits,basedir+'work/SPP/TestImage/WISPR-EM3_FM1_Inner.fits',hdr_Inner_0,img_Inner_0
  mreadfits,basedir+'work/SPP/TestImage/WISPR-EM3_FM1_Outer.fits',hdr_Outer_0,img_Outer_0

  suffix_binfac = ''
 ;Rebin by binfac if requested
  if keyword_set(bin) then begin
     if not keyword_set(bf) then bf=4
     suffix_binfac = '_binfac'+strmid(bf,7,1)
     img_Inner_0 = rebin(img_Inner_0,hdr_Inner_0.Naxis1/bf,hdr_Inner_0.Naxis2/bf)
     img_Outer_0 = rebin(img_Outer_0,hdr_Outer_0.Naxis1/bf,hdr_Outer_0.Naxis2/bf)
     hdr_Inner_0.Naxis1=hdr_Inner_0.Naxis1/bf
     hdr_Inner_0.Naxis2=hdr_Inner_0.Naxis2/bf
     hdr_Outer_0.Naxis1=hdr_Outer_0.Naxis1/bf
     hdr_Outer_0.Naxis2=hdr_Outer_0.Naxis2/bf
  endif

  suffix_square=''
 ;Create SQUARE images with the Test image located in its center,
 ;leaving all other pixels set to zero. Correct NAXIS2 in headers accordingly.
 ;This code uses the BINFACTOR corrected sizes above if /binfac was set.
  if keyword_set(squareFOV) then begin
     suffix_square='_squareFOV'
    ;Save original NAXIS2 in new variables.
     naxis2_inner_original = hdr_Inner_0.NAXIS2
     naxis2_outer_original = hdr_Inner_0.NAXIS2
    ;Correct NAXIS2 in headers to make square images.
     hdr_Inner_0.NAXIS2    = hdr_Inner_0.NAXIS1
     hdr_Outer_0.NAXIS2    = hdr_Outer_0.NAXIS1
    ;Create empty (-999) square images.
     tmp_Inner = fltarr(hdr_Inner_0.NAXIS1,hdr_Inner_0.NAXIS2) - 999.
     tmp_Outer = fltarr(hdr_Outer_0.NAXIS1,hdr_Outer_0.NAXIS2) - 999.
    ;Fit non-square images into center of (-999) square images.
     index0_inner = (hdr_Inner_0.NAXIS2-naxis2_inner_original)/2
     index0_outer = (hdr_Outer_0.NAXIS2-naxis2_outer_original)/2
     tmp_Inner(*,index0_inner:index0_inner+naxis2_inner_original-1) = img_Inner_0
     tmp_Outer(*,index0_outer:index0_outer+naxis2_outer_original-1) = img_Outer_0
     img_Inner_0 = tmp_Inner
     img_Outer_0 = tmp_Outer
  endif

  CircularOrbit_flag = 0
  experiment_suffix   = ''
  if keyword_set(CircularOrbits) then begin
     experiment_suffix='.CircularOrbits3degUnifStep_'
     CircularOrbit_flag = 1
  endif

 if keyword_set(UniformLong) then experiment_suffix = '.UniformLong'
  
  if keyword_set(FullList)  then begin
  listtype = 'full'
  filename = 'table_spp_orbits_full'+suffix_square+suffix_binfac+experiment_suffix+'.dat'
  endif

  if keyword_set(ShortList) then begin
  listtype = 'short'
  filename = 'table_spp_orbits_short'+suffix_square+suffix_binfac+experiment_suffix+'.dat'
  endif

  filename2  = 'table_relative_difference_rad_lon_lat'+suffix_square+suffix_binfac+experiment_suffix+'.dat'

 ;Set to -1. any pixel with value zero
  p = where(img_Inner_0 eq 0.) & if p[0] ne -1 then img_Inner_0(p)=-1.
  p = where(img_Outer_0 eq 0.) & if p[0] ne -1 then img_Outer_0(p)=-1.
 
 ;Store FOVs width and height [px] in convenient variables:
  FOVW_inner_px = hdr_Inner_0.NAXIS1
  FOVH_inner_px = hdr_Inner_0.NAXIS2
  FOVW_outer_px = hdr_Outer_0.NAXIS1
  FOVH_outer_px = hdr_Outer_0.NAXIS2

  ; set correction parameter
  if not keyword_set(correction) then abcorr = 'NONE'

  ; set orbital ephemerides:
  SciOrbNum   = [1,12,24]
  Norbits     = n_elements(SciOrbNum)

if keyword_set(SciOrbBrief) then begin
  Nepochs     = 3
  EpochArray  = strarr(Norbits,Nepochs) 
  ; These are approximate numbers, specially PERIHELION, which I took visually out of vizzer, refine using sobrbet,1,'SPP20180731P2_ephem.bsp'
  ;                  START                   PERIHELION              END
  EpochArray(0,*)  = ['2018-10-26  21:42:28','2018-11 01  03:42:28','2018-11-07  03:42:28']
  EpochArray(1,*)  = ['2022-05-23  04:47:11','2022-05-28  06:47:11','2022-06-02  11:47:11']
  EpochArray(2,*)  = ['2025-06-09  19:20:06','2025-06-14  16:20:06','2025-06-19  13:20:06']
; Convert the Epoch to ephemeris time (ET).
  ETarray = dblarr(Norbits,Nepochs)
  ETval=0.d
  for i=0,Norbits-1 do begin
  for j=0,Nepochs-1 do begin
     cspice_str2et, EpochArray[i,j], ETval
     ETarray[i,j] = ETval
  endfor
  endfor
endif

if keyword_set(ExtendedOrbits) or keyword_set(CircularOrbits) then begin
   if keyword_set(ExtendedOrbits) AND NOT keyword_set(ConstantCadence) AND NOT keyword_set(UniformLong) then begin
      print,'When selecting /ExtendeOrbits please specify: /ConstantCadence or /UniformLong.'
      return
   endif

   StartEpoch = ['2018-10-15  10:00:00','2022-05-12  10:00:00','2025-05-30  10:00:00']
   cspice_str2et, StartEpoch, StartET

if keyword_set(ConstantCadence) then begin
   NumDays         = 30.                          ; days 
   if NOT keyword_set(Cadence) then Cadence = 0.5 ; days
   CadenceSeconds  = Cadence * cspice_spd()       ; secs. It is the step in ET, in double precision. 
   Nepochs         = fix(Numdays / Cadence)
   EpochArray      = strarr(Norbits,Nepochs)
   ETarray         = dblarr(Norbits,Nepochs)
   ETarray   [*,0] = StartET[0:Norbits-1]
   for i=0,Norbits-1 do begin
       ETarray   [i,*] = ETarray[i,0] + CadenceSeconds * findgen(Nepochs)
       cspice_timout, reform(ETarray[i,*]), 'YYYY-MM-DD  HR:MN:SC', 22, EpochVector
       EpochArray[i,*] = EpochVector
    endfor
endif

if keyword_set(UniformLong) then begin
   CadenceSeconds0 = 6.0 * 3600.d0  ; sec  (6 hr cadence as base, it will be decreased when needed).
   Nepochs         = 1000           ; Excedengly large number 
   EpochArray      = strarr(Norbits,Nepochs)
   ETarray         = dblarr(Norbits,Nepochs) - 666.
   Longarray       = fltarr(Norbits,Nepochs) - 666.
   EpochArray[*,0] = StartEpoch
   ETarray   [*,0] = StartET
   for io=0,Norbits-1 do begin
      ET1 = reform(ETarray[io,0])
      cspice_timout, ET1, 'YYYY-MM-DD  HR:MN:SC', 22, EPOCH1
      DATE1   = strmid(EPOCH1,0,10)+'T'+strmid(EPOCH1,12,10)
      POS1    = get_sunspice_lonlat( DATE1, 'SPP', system='Carrington',/meters,/degrees)
      Longarray[io,0] = POS1[1] ; deg
   endfor
  ; Increace ET1 every CadenceSeconds and store ET1 and its EPOCH1 when
  ; achiveing a PSP sub-longitude step of DeltaLong within EPS accuracy
   EPS = 1.e-2
   if not keyword_set(DeltaLong) then DeltaLong = 3. ; deg
   for io=0,Norbits-1 do begin
      it  = 1
      continue_in_orbit:
      ET1            = reform(ETarray[io,it-1])
      CadenceSeconds = CadenceSeconds0
      add_CadenceSeconds:
      ET1 = ET1 + CadenceSeconds
      cspice_timout, ET1, 'YYYY-MM-DD  HR:MN:SC', 22, EPOCH1
      DATE1    = strmid(EPOCH1,0,10)+'T'+strmid(EPOCH1,12,10)
      POS1     = get_sunspice_lonlat( DATE1, 'SPP', system='Carrington',/meters,/degrees)
      SUBLON1  = POS1[1] ; deg
      DIST1    = POS1[0] ; m
      long1    = reform(Longarray[io,it-1])
      long2    = SUBLON1
      LongStep = long2 - long1
     ; Next line deals with the crossing of point Long=0.
      if long1 lt DeltaLong AND long2 gt 360.-DeltaLong then LongStep = LongStep - 360.
      if long2 lt DeltaLong AND long1 gt 360.-DeltaLong then LongStep = LongStep + 360.
     ;if LongStep lt 0. then goto,add_CadenceSeconds; Proceed again.
      if ABS(LongStep) lt DeltaLong*(1.-EPS) then goto,add_CadenceSeconds
      if ABS(LongStep) gt DeltaLong*(1.+EPS) then begin
         ET1 = ET1 - CadenceSeconds           ; Return to previous ET1 value.
         CadenceSeconds = CadenceSeconds / 2. ; Adapt CadenceSeconds if it was too large.
        ;print,'LongStep    ='+string(LongStep)      +' deg.'
        ;print,'New Cadence ='+string(CadenceSeconds)+' sec.'
         goto,add_CadenceSeconds              ; Proceed again.
      endif
      EpochArray[io,it] = EPOCH1
      ETarray   [io,it] = ET1
      Longarray [io,it] = SUBLON1
      it = it + 1
      print,'Orbit '+string(io)+string(it)+'   '+EPOCH1+string(ET1)+string(SUBLON1)+string(LongStep)
      if DIST1 le 0.51*AU then goto,continue_in_orbit
   endfor; Endorbit.
endif ; If UniformLong was set

if CircularOrbit_flag eq 1 then begin
   Norbits=1
   Radius = [10.,40.,80.]
   Nepochs = 120
   delta_carlon = 360./Nepochs  ; deg
   Cadence = 0.25 ; days
   CadenceSeconds  = Cadence * cspice_spd()       ; secs. It is the step in ET, in double precision. 
   EpochArray      = strarr(Norbits,Nepochs)
   ETarray         = dblarr(Norbits,Nepochs)
   ETarray   [*,0] = StartET[0:Norbits-1]
   for i=0,Norbits-1 do begin
       ETarray   [i,*] = ETarray[i,0] + CadenceSeconds * findgen(Nepochs)
       cspice_timout, reform(ETarray[i,*]), 'YYYY-MM-DD  HR:MN:SC', 22, EpochVector
       EpochArray[i,*] = EpochVector
    endfor
  
         if NOT keyword_set(Equatorial) AND NOT keyword_set(OffEquator) then begin
            print,'When selecting /CircularOrbits please specify: /Equatorial or /OffEquator.'
            return
         endif
         if keyword_set(Equatorial) then latitudes = fltarr(Nepochs)    
         if keyword_set(OffEquator) then begin
            amplitude_carlat = 3.           ; deg
            latitudes = -amplitude_carlat + 2.*amplitude_carlat*findgen(Nepochs)/float(Nepochs-1) ; deg
         endif        

   endif
endif

  science_limit       = 0.25      ; AU
  if CircularOrbit_flag eq 1 then science_limit = 1.0 ; AU
  flag_science_region = 0
; Make output table of ephemeris for the selected ETs
  openfile,1,filename
  openfile,2,filename2
  writekeys2,fileID=2
  count=1
  for i=0,Norbits-1 do begin
  writekeys,fileID=1
  if CircularOrbit_flag eq 1 then CircularOrbit_Radius = Radius[i]
  if keyword_set(UniformLong) then begin
     ET_vector = reform(ETarray(i,*))
     p = where(ET_vector ne -666.)
     Nepochs = n_elements(p)
  endif
  for j=0,Nepochs-1 do begin
  epoch = EpochArray[i,j]
     et = ETarray   [i,j]
  if CircularOrbit_flag eq 1 then begin
     CircularOrbit_Carlon = j*delta_carlon
     CircularOrbit_Carlat = latitudes[j]
  endif
  get_SPP_ephemeris,epoch=epoch,et=et,abcorr=abcorr ;,/printout
  print,'Done with date',count,' of',Norbits*Nepochs
  count=count+1
  if flag_science_region eq 0 AND dist_SUN_SPP/au le science_limit then begin
     flag_science_region = 1
     dot_line,fileID=1
  endif
  if flag_science_region eq 1 AND dist_SUN_SPP/au gt science_limit then begin
     flag_science_region = 0
     dot_line,fileID=1
  endif
  writedata,fileID=1,SciOrbNum=SciOrbNum,epocharray=epocharray,etarray=etarray
  if keyword_set(CreateFITS) then Create_FITS
  endfor; Loop Epochs of given Orbit
  endfor; Loop Orbits
  
  dash_line,fileID=1
  closefiles
 
  ; It's always good form to unload kernels after use, particularly in IDL due to data persistence.
  cspice_kclear

print,'Done.'

return
end

;   ephemeris_wispr,epoch='2018-10-26  21:42:28'
pro ephemeris_wispr,epoch=epoch,correction=correction,loadk=loadk,SciOrbNum=SciOrbNumb,FullList=FullList,ShortList=ShortList

common constants,c,rsun,au
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,px_inner_arcsec,px_outer_arcsec,sun_spp_vector_HCI,sun_spp_vector_HAE,sun_spp_vector_HEE,sun_spp_vector_HEQ
common output,listtype
common orbitnum,SciOrbNum

if keyword_set(FullList)  then begin
  listtype = 'full'
  filename = 'table_spp_orbits_full.dat'
endif

if keyword_set(ShortList) then begin
  listtype = 'short'
  filename = 'table_spp_orbits_short.dat'
endif

;load useful constants
 loadconstants

; LOAD KERNELS
  if keyword_set(loadk) then begin
      ;; Load kernels to initialize program.
      ;; Update by W. Thompson March 2017 to leverage new sunspice.
      load_sunspice_gen  
      ;; Load a set of kernels:  Use a meta kernel for convenience.
      cspice_furnsh, 'wispr_albert.tm'
  endif

  ; set correction parameter
  if not keyword_set(correction) then abcorr    = 'NONE'
  if not keyword_set(SciOrbNum)  then SciOrbNum = 0
  cspice_str2et, epoch, et
  get_SPP_ephemeris,epoch=epoch,et=et,abcorr=abcorr;,/printout
  writekeys,/terminal
  writedata,/terminal,SciOrbNum=SciOrbNum,epocharray=epoch,etarray=et
  dash_line,/terminal
  print
return
end  

pro get_SPP_ephemeris,epoch=epoch,et=et,abcorr=abcorr,printout=printout

common constants,c,rsun,au
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,px_inner_arcsec,px_outer_arcsec,sun_spp_vector_HCI,sun_spp_vector_HAE,sun_spp_vector_HEE,sun_spp_vector_HEQ
common FOVs,FOVW_inner_px,FOVH_inner_px,FOVW_outer_px,FOVH_outer_px
common experiment,CircularOrbit_flag,CircularOrbit_Radius,CircularOrbit_Carlon,CircularOrbit_Carlat,experiment_suffix
common sun_obs_unit_vector,sun_obs_unit_J2k,sun_obs_unit_HAE
common FOV_POINTINGS,alphaI_C,deltaI,alphaI_E,alphaI_W,alphaO_C,deltaO,alphaO_E,alphaO_W ; all in deg

    ;; Define parameters for a position lookup:
     ; Return the position vector of SPP as seen from SUN in the J2000 at the EPOCH
      observer = 'SUN'
      frame    = 'J2000'
;     frame    = 'HCI'

goto,skip_earth
      ;; Look-up the position for EARTH:
      target   = 'EARTH'
      cspice_spkpos, target, et, frame, abcorr, observer, ptarg, ltime
      ;; Output...
      print, 'The position of: '+target+', as observed from: '+observer+', in the reference frame: '+frame+', at EPOCH: '+EPOCH+', is:'
      print, '(x,y,z) (kilometers)     : ', ptarg[0:2]       ;; The first three entries of 'ptarg' contain the X, Y, Z position components.
      dist_SUN_SPP = sqrt(total(ptarg[0:2]^2))
      print, 'Distance Sun-'+target+': '+strtrim(dist_SUN_SPP/rsun,2)+' Rs = '+strtrim(dist_SUN_SPP,2)+' km = '+strtrim(dist_SUN_SPP/au,2)+' AU'
      print, 'Light ray length: '+strtrim(     c*ltime/rsun,2)+' Rs = '+strtrim(     c*ltime,2)+' km = '+strtrim(     c*ltime/au,2)+' AU'
      print
skip_earth:

      ;; Look-up the position for SPP:
      target   = 'SPP'
      cspice_spkpos, target, et, frame, abcorr, observer, ptarg, ltime
      sun_spp_vector_J2000 = ptarg[0:2] * 1.e3           ; m
      dist_SUN_SPP = sqrt(total(sun_spp_vector_J2000^2)) ; m
      ; cspice_vnorm(sun_spp_vector_J2000) gives the exact same result.

      spacecraft         = 'PSP';-96
      date               = strmid(epoch,0,10)+'T'+strmid(epoch,12,10)
      sun_spp_vector_HCI = get_sunspice_coord( date, spacecraft, system='HCI',/novelocity,/meters) ; m
      sun_spp_vector_HAE = get_sunspice_coord( date, spacecraft, system='HAE',/novelocity,/meters) ; m
      sun_spp_vector_HEE = get_sunspice_coord( date, spacecraft, system='HEE',/novelocity,/meters) ; m
      sun_spp_vector_HEQ = get_sunspice_coord( date, spacecraft, system='HEQ',/novelocity,/meters) ; m

if CircularOrbit_flag eq 1 then begin
   circular_experiment_position
   dist_SUN_SPP = CircularOrbit_Radius * rsun                                                      ; m
   sun_spp_vector_J2000 = sun_obs_unit_J2k * CircularOrbit_Radius * rsun                           ; m
   sun_spp_vector_HAE   = sun_obs_unit_HAE * CircularOrbit_Radius * rsun                           ; m
   ; Incorrectly set HCI, HEE, and HEQ, to J2k. They are never used anyhow:
   sun_spp_vector_HCI = sun_spp_vector_J2000
   sun_spp_vector_HEE = sun_spp_vector_J2000
   sun_spp_vector_HEQ = sun_spp_vector_J2000
endif 

      dist_SUN_SPP_HCI   = sqrt(total(sun_spp_vector_HCI^2))                                       ; m
      dist_SUN_SPP_HAE   = sqrt(total(sun_spp_vector_HAE^2))                                       ; m
      dist_SUN_SPP_HEE   = sqrt(total(sun_spp_vector_HEE^2))                                       ; m
      dist_SUN_SPP_HEQ   = sqrt(total(sun_spp_vector_HEQ^2))                                       ; m

;     Make sure all distances match within EPS accuracy
      EPS=1.e-5
      if abs(1.-dist_SUN_SPP/dist_SUN_SPP_HCI) gt EPS OR abs(1.-dist_SUN_SPP/dist_SUN_SPP_HAE) gt EPS OR abs(1.-dist_SUN_SPP/dist_SUN_SPP_HEE) gt EPS OR abs(1.-dist_SUN_SPP/dist_SUN_SPP_HEQ) gt EPS then begin
       print,'Coordinate results from cspice_spkpos and get_sunspice_coord do not match.'
       STOP
      endif
      
      if keyword_set(printout) then begin
      print, 'The position of: '+target+', as observed from: '+observer+', in the reference frame: '+frame+', at EPOCH: '+EPOCH+', is:'
      print, '(x,y,z) (kilometers)     : ', sun_spp_vector_J2000    ;; The first three entries of 'ptarg' contain the X, Y, Z position components.
      print, 'Distance Sun-'+target+': '+strtrim(dist_SUN_SPP/rsun,2)+' Rs = '+strtrim(dist_SUN_SPP,2)+' km = '+strtrim(dist_SUN_SPP/au,2)+' AU'
      print, 'Light ray length: '+strtrim(     c*ltime/rsun,2)+' Rs = '+strtrim(     c*ltime,2)+' km = '+strtrim(     c*ltime/au,2)+' AU'
      print
      endif

    ;Sub-spacecraft Carrington calculation
     target_ID = 'SPP'
 ;   target_ID = 'EARTH'
     frame = 'IAU_SUN'
     cspice_spkezr, target_ID, et, frame, abcorr, observer, kernel_targ1, kernel_ltime
						xcomp_au_start = kernel_targ1[0,0] / au
						ycomp_au_start = kernel_targ1[1,0] / au
						zcomp_au_start = kernel_targ1[2,0] / au
						Dist_xy_start=cspice_vnorm([xcomp_au_start, ycomp_au_start, 0])
						long_start = atan(ycomp_au_start,xcomp_au_start)*180/!Pi
						lat_start = atan(zcomp_au_start,Dist_xy_start)*180/!Pi
						vec = [0,0,0] ; Sun
                                                dis = get_dis( kernel_targ1[0:2,0], vec ) * 1.e3 ; m

                                                if long_start lt 0. then long_start = long_start + 360.
                                                
;  Recompute sub-SPP lat lon using Bill's routine
   frame = 'Carrington'
   r_CARR = get_sunspice_lonlat( date, spacecraft, system=frame,/meters,/degrees)
   
;  Make sure both results match within EPS accuracy
   printf,2,abs(1.-r_carr[0]/dis), abs(1.-r_carr[1]/long_start), abs(1.-r_carr[2]/lat_start)
   EPS=1.e-2
   if abs(1.-r_carr[0]/dis) gt EPS OR abs(1.-r_carr[1]/long_start) gt EPS OR abs(1.-r_carr[2]/lat_start) gt EPS then begin
      print,'Coordinate results from cspice_spkezr and get_sunspice_lonlat do not match.'
      STOP
   endif
 
;  Use Bill's routine output:
   dis        = r_CARR[0] ; m
   long_start = r_CARR[1] ; deg
   lat_start  = r_CARR[2] ; deg

   if CircularOrbit_flag eq 1 then begin
      dis        = CircularOrbit_Radius * rsun ; m
      long_start = CircularOrbit_Carlon        ; deg
      lat_start  = CircularOrbit_Carlat        ; deg
   endif
   
   						lon_string = STRTRIM(long_start,2)
						lat_string = STRTRIM(lat_start,2)
						IF long_start GE 0 then lon_string = '+' + lon_string
                                                IF lat_start  GE 0 then lat_string = '+' + lat_string
                                                
   goto,skip_test
   epoch_test = '2009-03-16  12:01:00.005'  
   cspice_str2et, epoch_test, et_test
   date_test = strmid(epoch_test,0,10)+'T'+strmid(epoch_test,12,13)
   spacecraft_test = 'STEREO_A'
   frame_test = 'Carr'
   r_CARR = get_sunspice_lonlat( date_test, spacecraft_test, system=frame_test,/meters)
   print,r_CARR[0],r_CARR[1:2]/!dtor
   STOP
   skip_test:
   
 if keyword_set(printout) then begin
   print,'Calculation of distance between '+OBSERVER+' and '+TARGET+', using frame '+frame+': '+strtrim(dis,2)+' km.'
   print,'Carrington Longitude: ' + lon_string + ' deg'
   print,'Latitude:             ' + lat_string + ' deg'
   print
 endif

      goto,skip
      SPP_WISPR_INNER_NAIF_ID = -96100
      room1                   =   4
      cspice_getfov, SPP_WISPR_INNER_NAIF_ID, room1, shape1, frame1, bsight1, bounds1

      SPP_WISPR_OUTER_NAIF_ID = -96120
      room2                   =   4
      cspice_getfov, SPP_WISPR_OUTER_NAIF_ID, room2, shape2, frame2, bsight2, bounds2
      skip:

      ; WISPRI/O swing (Z) and fov angles, LOADED 
      wispr_inner_Z_angle   = alphaI_C * !dtor ; rad
      wispr_inner_FOV_angle = deltaI   * !dtor ; rad
      wispr_outer_Z_angle   = alphaO_C * !dtor ; rad
      wispr_outer_FOV_angle = deltaO   * !dtor ; rad

      ; WISPRI/O angles for FOV East edge (E) and FOV West edge (W)
      wispr_inner_E_angle = wispr_inner_Z_angle - wispr_inner_FOV_angle / 2. ; rad
      wispr_inner_W_angle = wispr_inner_Z_angle + wispr_inner_FOV_angle / 2. ; rad

      wispr_outer_E_angle = wispr_outer_Z_angle - wispr_outer_FOV_angle / 2. ; rad
      wispr_outer_W_angle = wispr_outer_Z_angle + wispr_outer_FOV_angle / 2. ; rad

      ; PLANE-OF-SKY distance to the Sun Center of East, Center and West each FOV.
      D_SunCenter_inner_Center   = dist_SUN_SPP * sin(wispr_inner_Z_angle) ; m
      D_SunCenter_inner_EastEdge = dist_SUN_SPP * sin(wispr_inner_E_angle) ; m
      D_SunCenter_inner_WestEdge = dist_SUN_SPP * sin(wispr_inner_W_angle) ; m

      D_SunCenter_outer_Center   = dist_SUN_SPP * sin(wispr_outer_Z_angle) ; m
      D_SunCenter_outer_EastEdge = dist_SUN_SPP * sin(wispr_outer_E_angle) ; m
      D_SunCenter_outer_WestEdge = dist_SUN_SPP * sin(wispr_outer_W_angle) ; m

      ; Pixel size in arcsec:
      px_inner_arcsec = 3600.* (wispr_inner_FOV_angle/!dtor) / FOVW_inner_px
      px_outer_arcsec = 3600.* (wispr_outer_FOV_angle/!dtor) / FOVW_outer_px

      ; Pixel location of SunCenter in each FOV, assuming (ix,iy)=(0,0) in lower-left corner of each FOV.
      Pos_SunCenter_px_inner = [FOVW_inner_px / 2. + 1./2 - (wispr_inner_Z_angle/wispr_inner_FOV_angle) * FOVW_inner_px , FOVH_inner_px / 2 + 1./2]
      Pos_SunCenter_px_outer = [FOVW_outer_px / 2. + 1./2 - (wispr_outer_Z_angle/wispr_outer_FOV_angle) * FOVW_outer_px , FOVH_outer_px / 2 + 1./2]

if keyword_set(printout) then begin
      print,'--------------------------------------------------------------------'
      print,'Distance [Rsun] between Sun center and:'
      print,'--------------------------------------------------------------------'
      print,'                       East-Edge       Center          West-Edge'
      print,'WISPR INNER FOV:',[D_SunCenter_inner_EastEdge,D_SunCenter_inner_Center,D_SunCenter_inner_WestEdge]/rsun
      print,'WISPR OUTER FOV:',[D_SunCenter_outer_EastEdge,D_SunCenter_outer_Center,D_SunCenter_outer_WestEdge]/rsun
      print
      print,'--------------------------------------------------------------------'
      print,'Coordinates of Sun Center in pixels in each FOV, taking (0,0) in the lower left edge of each FOV:'
      print,'Inner:',Pos_SunCenter_px_inner
      print,'Outer:',Pos_SunCenter_px_outer
endif

Distances_SUN_FOV_inner_Rsun = [D_SunCenter_inner_EastEdge,D_SunCenter_inner_Center,D_SunCenter_inner_WestEdge]/rsun
Distances_SUN_FOV_outer_Rsun = [D_SunCenter_outer_EastEdge,D_SunCenter_outer_Center,D_SunCenter_outer_WestEdge]/rsun

return
end

pro openfile,n,filename
openw,n,filename
return
end

pro closefiles
close,/all
return
end

pro writedata,fileID=fileID,terminal=terminal,SciOrbNum=SciOrbNum,epocharray=epocharray,etarray=etarray
common output,listtype
common constants,c,rsun,au
common indexes,i,j
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,px_inner_arcsec,px_outer_arcsec,sun_spp_vector_HCI,sun_spp_vector_HAE,sun_spp_vector_HEE,sun_spp_vector_HEQ

if NOT keyword_set(terminal) then begin

if listtype eq 'full' then $  
  printf,fileID,SciOrbNum[i],epocharray[i,j],etarray[i,j],$
           sun_spp_vector_J2000/1.e3, dist_SUN_SPP/1.e3, dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Pos_SunCenter_px_inner, Pos_SunCenter_px_outer,$
           Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,$
           format='(1(I4),1(A22),1(E19.10),4(F14.2),2(F10.4),2(F10.2),4(F10.1),6(F10.4))'
if listtype eq 'short' then $  
  printf,fileID,SciOrbNum[i],epocharray[i,j],$
           dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Distances_SUN_FOV_inner_Rsun([0,2]), Distances_SUN_FOV_outer_Rsun([0,2]),$
           format='(1(I4),1(A22),2(F10.4),2(F10.2),4(F10.4))'

endif else begin

if listtype eq 'full' then $  
  print   ,SciOrbNum,epocharray,etarray,$
           sun_spp_vector_J2000/1.e3, dist_SUN_SPP/1.e3, dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Pos_SunCenter_px_inner, Pos_SunCenter_px_outer,$
           Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,$
           format='(1(I4),1(A22),1(E19.10),4(F14.2),2(F10.4),2(F10.2),4(F10.1),6(F10.4))'
if listtype eq 'short' then $  
  print   ,SciOrbNum,epocharray,$
           dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Distances_SUN_FOV_inner_Rsun([0,2]), Distances_SUN_FOV_outer_Rsun([0,2]),$
           format='(1(I4),1(A22),2(F10.4),2(F10.2),4(F10.4))'
endelse

return
end

pro writekeys,fileID=fileID,terminal=terminal
common output,listtype

if listtype eq 'full'  then begin  
string1=' Orb     DATE       TIME           ET            SObsJ2000x    SObsJ2000y    SObsJ2000z        D-SObs    D-SObs    D-SObs       LON       LAT    CtrI_x    CtrI_y    CtrO_x    CtrO_y      S-IE      S-IC      S-IW      S-OE      S-OC      S-OW'
string2='   #                                                   [km]          [km]          [km]          [km]    [Rsun]      [AU]     [deg]     [deg]      [px]      [px]      [px]      [px]    [Rsun]    [Rsun]    [Rsun]    [Rsun]    [Rsun]    [Rsun]' 
endif else begin
string1=' Orb     DATE       TIME      D-SObs    D-SObs       LON       LAT      S-IE      S-IW      S-OE      S-OW'
string2='   #                          [Rsun]      [AU]     [deg]     [deg]    [Rsun]    [Rsun]    [Rsun]    [Rsun]' 
endelse

if not keyword_set(terminal) then begin
dash_line,fileID=fileID
printf,fileID,string1
printf,fileID,string2
dash_line,fileID=fileID
endif else begin
dash_line,/terminal
print,string1
print,string2
dash_line,/terminal
endelse

return
end

pro writekeys2,fileID=fileID
printf,fileID,'FracDiff(r),FracDiff(Lon),FracDiff(Lat)'
return
end

pro dash_line,fileID=fileID,terminal=terminal
  common output,listtype
if listtype eq 'full'  then $
string1='--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
if listtype eq 'short' then $
string1='----------------------------------------------------------------------------------------------------------'
if not keyword_set(terminal) then printf,fileID,string1
if     keyword_set(terminal) then print        ,string1
return
end

pro dot_line,fileID=fileID
common output,listtype
if listtype eq 'full'  then $
printf,fileID,'  ................................................................................................................................................................................................................................................'
if listtype eq 'short' then $
printf,fileID,'  .......................................................................................................'
return
end

pro loadconstants
common constants,c,rsun,au
; Set useful constants
  c    = 299792.458e3    ; m/sec
  rsun = 695700.e3       ; m 
  au   = 149.597870700e9 ; m
return
end

FUNCTION get_dis, camera_xyz, vec
	dis = SQRT( (camera_xyz[0]-vec[0])^2 + (camera_xyz[1]-vec[1])^2 + (camera_xyz[2]-vec[2])^2 )
	RETURN,dis
END

pro Create_FITS
common constants,c,rsun,au
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,px_inner_arcsec,px_outer_arcsec,sun_spp_vector_HCI,sun_spp_vector_HAE,sun_spp_vector_HEE,sun_spp_vector_HEQ
common SynthFITS,hdr_Inner_0,hdr_Outer_0,img_Inner_0,img_Outer_0,datadir,epoch,et
common suffixes,suffix_binfac,suffix_square
common experiment,CircularOrbit_flag,CircularOrbit_Radius,CircularOrbit_Carlon,CircularOrbit_Carlat,experiment_suffix

; Create HDR equal to "TestImages" Header
hdr_Inner = hdr_Inner_0
hdr_Outer = hdr_Outer_0

; Update .DATE_OBS
hdr_Inner.DATE_OBS =  strmid(epoch,0,10)+'T'+strmid(epoch,12,10)
hdr_Outer.DATE_OBS =  strmid(epoch,0,10)+'T'+strmid(epoch,12,10)

; Create filenames
filename_Inner = 'WISPR_I_'+hdr_Inner.DATE_OBS+suffix_square+suffix_binfac+'_Blank'+experiment_suffix+'.fts'
filename_Outer = 'WISPR_O_'+hdr_Inner.DATE_OBS+suffix_square+suffix_binfac+'_Blank'+experiment_suffix+'.fts'

; Update .FILENAME
hdr_Inner.FILENAME = filename_Inner
hdr_Outer.FILENAME = filename_Outer

; Update .DSUN_OBS, distance [m] to Sun Center
hdr_Inner.DSUN_OBS = dist_SUN_SPP ; m
hdr_Outer.DSUN_OBS = dist_SUN_SPP ; m

; Update .CRPIX1
hdr_Inner.CRPIX1 = Pos_SunCenter_px_inner[0]
hdr_Outer.CRPIX1 = Pos_SunCenter_px_outer[0]

; Update .CRPIX2
hdr_Inner.CRPIX2 = Pos_SunCenter_px_inner[1]; the original headers say 961.5, ask Angelos why?
hdr_Outer.CRPIX2 = Pos_SunCenter_px_outer[1]

; Update .CDELT1
hdr_Inner.cdelt1 = px_inner_arcsec
hdr_Outer.cdelt1 = px_outer_arcsec

; Expand HEADERS with other needed variables

hdr_Inner = create_struct(hdr_Inner,$
            'CRLN_OBS'  ,long_start,$
            'CRLT_OBS'  , lat_start,$
            'CROTA'     , 0.,$
    	    'J2kX_OBS'  ,sun_spp_vector_J2000[0],$
    	    'J2kY_OBS'  ,sun_spp_vector_J2000[1],$
            'J2kZ_OBS'  ,sun_spp_vector_J2000[2],$
            'HCIX_OBS'  ,sun_spp_vector_HCI  [0],$
            'HCIY_OBS'  ,sun_spp_vector_HCI  [1],$
            'HCIZ_OBS'  ,sun_spp_vector_HCI  [2],$
            'HAEX_OBS'  ,sun_spp_vector_HAE  [0],$
            'HAEY_OBS'  ,sun_spp_vector_HAE  [1],$
            'HAEZ_OBS'  ,sun_spp_vector_HAE  [2],$
            'HEEX_OBS'  ,sun_spp_vector_HEE  [0],$
            'HEEY_OBS'  ,sun_spp_vector_HEE  [1],$
            'HEEZ_OBS'  ,sun_spp_vector_HEE  [2],$
            'HEQX_OBS'  ,sun_spp_vector_HEQ  [0],$
            'HEQY_OBS'  ,sun_spp_vector_HEQ  [1],$
            'HEQZ_OBS'  ,sun_spp_vector_HEQ  [2] )

hdr_Outer = create_struct(hdr_Outer,$
            'CRLN_OBS'  ,long_start,$
            'CRLT_OBS'  , lat_start,$
            'CROTA'     , 0.,$
    	    'J2kX_OBS'  ,sun_spp_vector_J2000[0],$
    	    'J2kY_OBS'  ,sun_spp_vector_J2000[1],$
            'J2kZ_OBS'  ,sun_spp_vector_J2000[2],$
            'HCIX_OBS'  ,sun_spp_vector_HCI  [0],$
            'HCIY_OBS'  ,sun_spp_vector_HCI  [1],$
            'HCIZ_OBS'  ,sun_spp_vector_HCI  [2],$
            'HAEX_OBS'  ,sun_spp_vector_HAE  [0],$
            'HAEY_OBS'  ,sun_spp_vector_HAE  [1],$
            'HAEZ_OBS'  ,sun_spp_vector_HAE  [2],$
            'HEEX_OBS'  ,sun_spp_vector_HEE  [0],$
            'HEEY_OBS'  ,sun_spp_vector_HEE  [1],$
            'HEEZ_OBS'  ,sun_spp_vector_HEE  [2],$
            'HEQX_OBS'  ,sun_spp_vector_HEQ  [0],$
            'HEQY_OBS'  ,sun_spp_vector_HEQ  [1],$
            'HEQZ_OBS'  ,sun_spp_vector_HEQ  [2] )

; Write FITS files to disk
mwritefits,hdr_Inner,img_Inner_0,outfile=datadir+filename_Inner
mwritefits,hdr_Outer,img_Outer_0,outfile=datadir+filename_Outer

return
end

; radlon_coverage_plot,table_file='table.ScienceOrbit01.short.txt'
; radlon_coverage_plot,table_file='table.ScienceOrbit12.short.txt'
; radlon_coverage_plot,table_file='table.ScienceOrbit24.short.txt'

pro wrapper_radlon_coverage

 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.01.UPDATED-POINTINGS.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.12.UPDATED-POINTINGS.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.24.UPDATED-POINTINGS.txt',/fov_edge_lon,/xtit
 
 return

 radlon_coverage_plot,table_file='table.CircularOrbit01.short.UPDATED-POINTINGS.txt',/fov_edge_lon,/xtit
 return
 
 
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.12.UPDATED-POINTINGS.txt',/fov_edge_lon



 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.24.UPDATED-POINTINGS.txt',/fov_edge_lon
 return
 

 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.01.UPDATED-POINTINGS.txt',/fov_edge_lon


 return
 
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.01.txt',/sub_psp_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.01.txt',/sub_psp_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.12.txt',/sub_psp_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.12.txt',/sub_psp_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.24.txt',/sub_psp_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.24.txt',/sub_psp_lon

return
  
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.01.txt',/fov_center_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.01.txt',/fov_center_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.12.txt',/fov_center_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.12.txt',/fov_center_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.24.txt',/fov_center_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.24.txt',/fov_center_lon

 return
 
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.24.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.24.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.12.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.12.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.SciOrbit.01.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.ExtOrbit.01.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.CircularOrbit01.short.txt',/fov_edge_lon
 return

 radlon_coverage_plot,table_file='table.UnifLong.MidOrbit.01.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.MidOrbit.12.txt',/fov_edge_lon
 radlon_coverage_plot,table_file='table.UnifLong.MidOrbit.24.txt',/fov_edge_lon
 
 return
 
 return
end

pro radlon_coverage_plot,table_file=table_file,input_dir=input_dir,$
                         fov_center_lon=fov_center_lon,fov_edge_lon=fov_edge_lon,sub_psp_lon=sub_psp_lon,$
                         xtit=xtit
  common FOV_POINTINGS,alphaI_C,deltaI,alphaI_E,alphaI_W,alphaO_C,deltaO,alphaO_E,alphaO_W ; all in deg

  load_fov_pointings

  if not keyword_set(input_dir) then input_dir='/data1/work/SPP/SORBET_VIZZER_WISPR_RevA/'

  readcol,input_dir+table_file,orb,date,time,Dsun_Rsun,Dsun_AU,Lon,Lat,HIE,HIW,HOE,HOW,$
          format='U,A,A,F,F,F,F,F,F,F,F',skipline=4,/quick

  Ndat = (size(Lon))(1)

  p = where(lon gt 360.)
  if p(0) ne -1 then stop

  if Orb[0] lt 10 then Orb_string = '0'+strmid(string(Orb[0]),7,1)
  if Orb[0] ge 10 then Orb_string =     strmid(string(Orb[0]),6,2)

  if keyword_set(sub_psp_lon)    then begin
     suffix = '_sub-psp-lon'
     xtitle = 'Sub-PSP Longitude [deg]'
     xmin   = - 30. ; deg 
     xmax   = +460. ; deg 
  endif
  
  if keyword_set(fov_edge_lon)   then begin
     suffix = '_fov-edge-lon'
     xtitle=''
     if keyword_set(xtit) then xtitle = 'FOV-Edges Longitude [deg]'
     xmin   =    0.     ; deg 
     xmax   = +360.+90. ; deg 
 endif
  
  if keyword_set(fov_center_lon) then begin
     suffix = '_fov-center-lon'
     xtitle = 'FOV-Center Longitude [deg]'
     xmin   = -100. ; deg 
     xmax   = +500. ; deg 
  endif

  ps1,input_dir+table_file+'_coverage_plot'+suffix+'.eps',0
  DEVICE,XSIZE=20.0,YSIZE=10.0,scale_factor=10
  !p.charsize=1.5

  dextlon   = 60
  Nextlon   = (xmax-xmin)/dextlon + 1
   extlon   = xmin+dextlon*findgen(Nextlon)
  extlonlab = strarr(Nextlon)
  for i=0,Nextlon-1 do begin
     tmp = extlon(i)
     if tmp lt   0. then tmp = tmp + 360.
     if tmp gt 360. then tmp = tmp - 360.
     if tmp ge   0. and tmp lt  10. then extlonlab(i) = strmid(string(tmp),6,1)
     if tmp ge  10. and tmp lt 100. then extlonlab(i) = strmid(string(tmp),6,2)
     if tmp ge 100. and tmp le 360. then extlonlab(i) = strmid(string(tmp),6,3)
  endfor
  
  exthow=max(how)+fltarr(Nextlon)

  goto,skip
  plot,lon,HOW,xtitle=xtitle,ytitle='Heliocentric Height [R!DSUN!N]',$
       title='Orbit-'+Orb_string+' FOVs: Inner (blue), Outer (red).',$
       xr=[xmin,xmax],xstyle=1,$
       font=1,/nodata
  skip:

  plot,extlon,exthow,$
       xtitle=xtitle,ytitle='Heliocentric Height [R!DSUN!N]',$
;      title='Circular Orbit of Radius 10 R!DSUN!N',yr=[1.,11.],ystyle=1,$
       title='PSP Orbit '+Orb_string,$     ;+' FOVs: Inner (blue), Outer (red).',$
       xr=[xmin,xmax],xstyle=1,$
       font=1,/nodata,$
       xticks=Nextlon,xtickname = extlonlab,xtickv=extlon

  oplot,  0*[1,1],max(how)*[0,2]
  oplot,360*[1,1],max(how)*[0,2]
  
  loadct,12
  blue=100
  red =200
  thick=1
  size=1
  perihelion = median(where(HIE eq min(HIE)))
;common FOV_POINTINGS,alphaI_C,deltaI,alphaI_E,alphaI_W,alphaO_C,deltaO,alphaO_E,alphaO_W ; all in deg  
  for i=0,Ndat-1 do begin
     lnsty=2
     
    if keyword_set(fov_edge_lon) then begin
     shift_IO =  0. ; deg
     deltaI_E = 90. - alphaI_E  & LonI_E = Lon[i] + deltaI_E  
     deltaI_C = 90. - alphaI_C  & LonI_C = Lon[i] + deltaI_C
     deltaI_W = 90. - alphaI_W  & LonI_W = Lon[i] + deltaI_W     
     deltaO_E = 90. - alphaO_E  & LonO_E = Lon[i] + deltaO_E + shift_IO
     deltaO_C = 90. - alphaO_C  & LonO_C = Lon[i] + deltaO_C + shift_IO
     deltaO_W = 90. - alphaO_W  & LonO_W = Lon[i] + deltaO_W + shift_IO
   ; Correct LonO_W with the idea of highlighting the max height that is
   ; best constrained by the LOSs.
     LonO_W = Lon[i]
   ; But then LonO_C needs to be corrected accordingly:
     LonO_C = (LonO_E+LonO_W)/2.
   ; On same spirit, assign HOW = Dsun_Rsun.
     HOW[i] = Dsun_Rsun[i]
   
  endif

    if keyword_set(fov_center_lon) then begin
       shift_IO =  0. ; deg
       deltaI_C = 90. - alphaI_C  & LonI_C = Lon[i] + deltaI_C
       deltaO_C = 90. - alphaO_C  & LonO_C = Lon[i] + deltaO_C + shift_IO
       LonI_E = LonI_C
       LonI_C = LonI_C
       LonI_W = LonI_C
       LonO_E = LonO_C + shift_IO
       LonO_C = LonO_C + shift_IO       
       LonO_W = LonO_C + shift_IO
    endif

    if keyword_set(sub_psp_lon) then begin
       shift_IO = 0.5 ; deg
       LonI_E = Lon[i]
       LonI_C = Lon[i]
       LonI_W = Lon[i]
       LonO_E = Lon[i] + shift_IO
       LonO_C = Lon[i] + shift_IO       
       LonO_W = Lon[i] + shift_IO
    endif
   
;    if i eq 0 or i eq perihelion or i eq Ndat-1 then thick=4
    if  Orb[0] eq 1                                                    then lnsty=0
    if (Orb[0] eq 12 or Orb[0] eq 24) AND (i+1 ge 7 AND i+1 le Ndat-5) then lnsty=0
     oplot,[LonI_E,LonI_W],[HIE[i],HIW[i]],linestyle=lnsty,color=blue
     oplot,[LonO_E,LonO_W],[HOE[i],HOW[i]],linestyle=lnsty,color=red
     if i eq 0 then begin
     oplot,[LonI_E],[HIE[i]],th=thick,color=blue,psym=4,symsize=size
     oplot,[LonO_W],[HOW[i]],th=thick,color=red ,psym=4,symsize=size
     endif
     if i eq perihelion then begin
     oplot,[LonI_E],[HIE[i]],th=thick,color=blue,psym=2,symsize=size
     oplot,[LonO_W],[HOW[i]],th=thick,color=red ,psym=2,symsize=size
     endif
     if i eq Ndat-1 then begin
     oplot,[LonI_E],[HIE[i]],th=thick,color=blue,psym=5,symsize=size
     oplot,[LonO_W],[HOW[i]],th=thick,color=red ,psym=5,symsize=size
     endif
  endfor
  loadct,0
  ps2
return
end

pro circular_experiment_position
common sun_obs_unit_vector,sun_obs_unit_J2k,sun_obs_unit_HAE
common experiment,CircularOrbit_flag,CircularOrbit_Radius,CircularOrbit_Carlon,CircularOrbit_Carlat,experiment_suffix

; /* J2000 solar pole coords */
  ALPHApo = 286.13 * !dtor
  DELTApo =  63.87 * !dtor

; /* solar pole vector */
  spol1    = dblarr(3)
  spol1[0] = cos(DELTApo) * cos(ALPHApo)
  spol1[1] = cos(DELTApo) * sin(ALPHApo)
  spol1[2] = sin(DELTApo)

; g1 unit vector:  g1 = spol1 x unit_x
; being unit_x = (1,0,0) the unit x vector of the J2000 coordinate system.
  g1     =  dblarr(3)
  g1[0]  =  0.
  g1[1]  =  spol1[2]
  g1[2]  = -spol1[1]
  normg1 =  sqrt(total(g1^2))
  g1     =  g1 / normg1

;  /* the J2000.0 angle between the Ecliptic and mean Equatorial planes
;   * is 23 deg 26 min 21.4119 sec - From Allen's Astrophysical
;                                    Quantities, 4th ed. (2000) */
; The added Carlat must be set to ZERO to get an equatorial orbit
  
  alpha = (23. + 26./60. + 21.4119/3600. + CircularOrbit_Carlat)*!dtor
  
; g0, being g1 in HAE coordinates.
  g0    =  dblarr(3)
  g0[0] =  g1[0]
  g0[1] =  cos(alpha)*g1[1] + sin(alpha)*g1[2]
  g0[2] = -sin(alpha)*g1[1] + cos(alpha)*g1[2]

; Store the unit position vector in J2k and HAE
  sun_obs_unit_J2k = g1
  sun_obs_unit_HAE = g0

end

pro load_fov_pointings
  common FOV_POINTINGS,alphaI_C,deltaI,alphaI_E,alphaI_W,alphaO_C,deltaO,alphaO_E,alphaO_W ; all in deg

goto,updated_values
; These are the original values we used for the Simulations carried
; out at CLASP in 2017 and Buenos Aires in early 2018, updated now by
; the numbers below.

  alphaI_C = 32.2 ; deg 
  deltaI   = 40.9 ; deg 
  alphaO_C = 77.0 ; deg 
  deltaO   = 59.2 ; deg 
  
  alphaI_E = alphaI_C - deltaI/2.
  alphaI_W = alphaI_C + deltaI/2.
  alphaO_E = alphaO_C - deltaO/2.
  alphaO_W = alphaO_C + deltaO/2.
  
  updated_values:
; These are the values updated by Angelos on February-21-2018, during
; the WISPR planning meeting at Maryland.
  alphaI_C = 0.5*( 51.4 + 13.5) ; deg 
  deltaI   =     ( 51.4 - 13.5) ; deg 
  alphaO_C = 0.5*(104.7 + 49.7) ; deg 
  deltaO   =     (104.7 - 49.7) ; deg 
  
  alphaI_E = alphaI_C - deltaI/2.
  alphaI_W = alphaI_C + deltaI/2.
  alphaO_E = alphaO_C - deltaO/2.
  alphaO_W = alphaO_C + deltaO/2.
  
  return
end

pro compare_spp_pointings
  common FOV_POINTINGS,alphaI_C,deltaI,alphaI_E,alphaI_W,alphaO_C,deltaO,alphaO_E,alphaO_W ; all in deg

  load_fov_pointings
  load_updated_fov_pointings

  v  = [alphaI_E,alphaI_W,alphaO_E,alphaO_W]
  v1 = [alphaI_E_1,alphaI_W_1,alphaO_E_1,alphaO_W_1]
  print,'WISPR pointing information:'
  print,'      T1_inner     T1_outer     T2_inner     T2_outter'
  print,'Used:    ',v
  print,'Updated: ',v1
  print,'Fractional relative difference of Used values relative to updated ones:'
  frdiff = (v-v1) / v1
  print,'          ',frdiff
  print
  print,'Pxiel size [arcsec]'
  print,'Used:   ',[deltaI  ,deltaO  ]*3600./2048
  print,'Updated:',[deltaI_1,deltaO_1]*3600./2048
  print,'Fractional relative difference of Used values relative to updated ones:'
  frdiff = ([deltaI  ,deltaO] - [deltaI_1,deltaO_1])/[deltaI_1,deltaO_1]
  print,'          ',frdiff
  print

  stop
  return
end
