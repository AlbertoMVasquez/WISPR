;;
;
; IDL tool to generate WISPR orbital ephemeris with relevant parameters
; for tomography. It also generates synthetic headers for the
; corresponding fits files.
;
; INPUTS:  w.i.p.
;
; OUTPUTS: w.i.p.
;  
; CALL SEQUENCE: wispr_tool, /loadk
;
; HISTORY
; A.M. Vasquez -  Oct-15-2017 - Version 1.0
;
;;

pro wispr_tool,loadk=loadk,correction=correction,SciOrbBrief=SciOrbBrief,ExtendedOrbits=ExtendedOrbits,FullList=FullList,ShortList=ShortList,CreateFITS=CreateFITS,Outdir=Outdir,basedir=basedir

common constants,c,rsun,au
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun, Distances_SUN_FOV_inner_px, Distances_SUN_FOV_outer_px,px_inner_arcsec,px_outer_arcsec
common output,listtype
common indexes,i,j
common SynthFITS,hdr_Inner_0,hdr_Outer_0,img_Inner_0,img_Outer_0,datadir,epoch,et

if keyword_set(FullList)  then begin
  listtype = 'full'
  filename = 'table_spp_orbits_full.dat'
endif

if keyword_set(ShortList) then begin
  listtype = 'short'
  filename = 'table_spp_orbits_short.dat'
endif

if NOT keyword_set(Outdir) then Outdir='../TestImage/'

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

  ; Load sample FITS if CreateFITS
  if keyword_set(CreateFITS) then begin
     if not keyword_set(basedir) then basedir='/data1/'
     mreadfits,basedir+'work/SPP/TestImage/WISPR-EM3_FM1_Inner.fits',hdr_Inner_0,img_Inner_0
     mreadfits,basedir+'work/SPP/TestImage/WISPR-EM3_FM1_Outer.fits',hdr_Outer_0,img_Outer_0
     img_Inner_0 = 0. * img_Inner_0
     img_Outer_0 = 0. * img_Outer_0
     datadir=Outdir
  endif
     
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
    cspice_str2et, EpochArray[i,j], ETval & ETarray[i,j] = ETval
  endfor
  endfor
endif

if keyword_set(ExtendedOrbits) then begin
   StartEpoch = ['2018-10-15  10:00:00','2022-05-12  10:00:00','2025-05-30  10:00:00']
   cspice_str2et, StartEpoch, StartET
   NumDays         = 30.                    ; days 
   Cadence         = 0.5                    ; days
   CadenceSeconds  = Cadence * cspice_spd() ; secs. It is the step in ET, in double precision. 
   Nepochs         = fix(Numdays / Cadence)
   EpochArray      = strarr(Norbits,Nepochs)
   ETarray         = dblarr(Norbits,Nepochs)
   ETarray   [*,0] = StartET
   for i=0,Norbits-1 do begin
       ETarray   [i,*] = ETarray[i,0] + CadenceSeconds * findgen(Nepochs)
       cspice_timout, reform(ETarray[i,*]), 'YYYY-MM-DD  HR:MN:SC', 22, EpochVector
       EpochArray[i,*] = EpochVector
   endfor
endif

  science_limit       = 0.25 ; AU
  flag_science_region = 0 
; Make output table of ephemeris for the selected ETs
  openfile,filename
  for i=0,Norbits-1 do begin
  writekeys
  for j=0,Nepochs-1 do begin  
  epoch = EpochArray[i,j]
     et = ETarray   [i,j]
  get_SPP_ephemeris,epoch=epoch,et=et,abcorr=abcorr;,/printout
  if flag_science_region eq 0 AND dist_SUN_SPP/au le science_limit then begin
     flag_science_region = 1
     point_line
  endif
  if flag_science_region eq 1 AND dist_SUN_SPP/au gt science_limit then begin
     flag_science_region = 0
     point_line
  endif
  writedata,SciOrbNum=SciOrbNum,epocharray=epocharray,etarray=etarray
  if keyword_set(CreateFITS) then Create_FITS
  endfor
  endfor
  dash_line
  closefile
 
  ; It's always good form to unload kernels after use, particularly in IDL due to data persistence.
  cspice_kclear

print,'Done.'

return
end

;   ephemeris_wispr,epoch='2018-10-26  21:42:28'
pro ephemeris_wispr,epoch=epoch,correction=correction,loadk=loadk,SciOrbNum=SciOrbNumb,FullList=FullList,ShortList=ShortList

common constants,c,rsun,au
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun, Distances_SUN_FOV_inner_px, Distances_SUN_FOV_outer_px,px_inner_arcsec,px_outer_arcsec
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
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun, Distances_SUN_FOV_inner_px, Distances_SUN_FOV_outer_px,px_inner_arcsec,px_outer_arcsec

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
      sun_spp_vector_J2000 = ptarg[0:2] ; km
      dist_SUN_SPP = sqrt(total(sun_spp_vector_J2000^2)) ; cspice_vnorm(sun_spp_vector_J2000) gives the exact same result.

      spacecraft = -96
      date       = strmid(epoch,0,10)+'T'+strmid(epoch,12,10)
      frame      = 'HAE'
      state      = get_sunspice_coord( date, spacecraft, system=frame)
      
      if keyword_set(printout) then begin
      print, 'The position of: '+target+', as observed from: '+observer+', in the reference frame: '+frame+', at EPOCH: '+EPOCH+', is:'
      print, '(x,y,z) (kilometers)     : ', sun_spp_vector_J2000    ;; The first three entries of 'ptarg' contain the X, Y, Z position components.
      print, 'Distance Sun-'+target+': '+strtrim(dist_SUN_SPP/rsun,2)+' Rs = '+strtrim(dist_SUN_SPP,2)+' km = '+strtrim(dist_SUN_SPP/au,2)+' AU'
      print, 'Light ray length: '+strtrim(     c*ltime/rsun,2)+' Rs = '+strtrim(     c*ltime,2)+' km = '+strtrim(     c*ltime/au,2)+' AU'
      print
      endif

    ;Sub-spacecraft Carrington calculation
     target_ID = 'SPP'
;    target_ID = 'EARTH'
     frame = 'IAU_SUN'
     cspice_spkezr, target_ID, et, frame, abcorr, observer, kernel_targ1, kernel_ltime
						xcomp_au_start = kernel_targ1[0,0] / au
						ycomp_au_start = kernel_targ1[1,0] / au
						zcomp_au_start = kernel_targ1[2,0] / au
						Dist_xy_start=cspice_vnorm([xcomp_au_start, ycomp_au_start, 0])
						long_start = atan(ycomp_au_start,xcomp_au_start)*180/!Pi
						lat_start = atan(zcomp_au_start,Dist_xy_start)*180/!Pi
						lon_string = STRTRIM(long_start,2)
						lat_string = STRTRIM(lat_start,2)
						IF long_start GE 0 then lon_string = '+' + lon_string
						IF lat_start  GE 0 then lat_string = '+' + lat_string
						vec = [0,0,0] ; Sun
						dis = get_dis( kernel_targ1[0:2,0], vec )

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

      wispr_inner_Z_angle = 32.2 * !dtor ; rad
      wispr_outer_Z_angle = 77.0 * !dtor ; rad
      
      D_SunCenter_inner_Center = dist_SUN_SPP * sin(wispr_inner_Z_angle) ; km
      D_SunCenter_outer_Center = dist_SUN_SPP * sin(wispr_outer_Z_angle) ; km

      wispr_inner_FOV_angle = 40.9 * !dtor ; rad
      wispr_outer_FOV_angle = 59.2 * !dtor ; rad

      alpha_inner = wispr_inner_Z_angle
       beta_inner = alpha_inner - wispr_inner_FOV_angle / 2.
      gamma_inner = !pi/2.      + wispr_inner_FOV_angle / 2.
      D_SunCenter_inner_EastEdge = dist_SUN_SPP * sin(beta_inner) / sin(gamma_inner)
      D_SunCenter_inner_WestEdge = 2.* D_SunCenter_inner_Center - D_SunCenter_inner_EastEdge 

      alpha_outer = wispr_outer_Z_angle
       beta_outer = alpha_outer - wispr_outer_FOV_angle / 2.
      gamma_outer = !pi/2.      + wispr_outer_FOV_angle / 2.
      D_SunCenter_outer_EastEdge = dist_SUN_SPP * sin(beta_outer) / sin(gamma_outer)
      D_SunCenter_outer_WestEdge = 2.* D_SunCenter_outer_Center - D_SunCenter_outer_EastEdge 

; FOVs width and height in Pixels:
  FOVW_px = 2048.
  FOVH_px = 1920.
  
; FOVs widths in Rsun:
  FOVW_inner_rsun = (D_SunCenter_inner_WestEdge-D_SunCenter_inner_EastEdge)/rsun
  FOVW_outer_rsun = (D_SunCenter_outer_WestEdge-D_SunCenter_outer_EastEdge)/rsun

; Pixel size in arcsec:
  px_inner_arcsec = 3600.* (wispr_inner_FOV_angle/!dtor) / FOVW_px
  px_outer_arcsec = 3600.* (wispr_outer_FOV_angle/!dtor) / FOVW_px

; Pixel size in Rsun:  
  px_inner_rsun = FOVW_inner_rsun / FOVW_px
  px_outer_rsun = FOVW_outer_rsun / FOVW_px

; Distances SunCenter-FOVcenter in pixels:
  D_SunCenter_inner_Center_px = D_SunCenter_inner_Center / rsun / px_inner_rsun
  D_SunCenter_outer_Center_px = D_SunCenter_outer_Center / rsun / px_outer_rsun

; Distances SunCenter-FOVEastEdge in pixels:
  D_SunCenter_inner_EastEdge_px = D_SunCenter_inner_EastEdge /rsun / px_inner_rsun
  D_SunCenter_outer_EastEdge_px = D_SunCenter_outer_EastEdge /rsun / px_outer_rsun

; Distances SunCenter-FOVEastEdge in pixels:
  D_SunCenter_inner_WestEdge_px = D_SunCenter_inner_WestEdge /rsun / px_inner_rsun
  D_SunCenter_outer_WestEdge_px = D_SunCenter_outer_WestEdge /rsun / px_outer_rsun

; Pixel location of SunCenter in each FOV, assuming (ix,iy)=(0,0) in
; lower-left corner of each FOV.
  Pos_SunCenter_px_inner = [- D_SunCenter_inner_EastEdge_px , FOVH_px / 2 + 1./2]
  Pos_SunCenter_px_outer = [- D_SunCenter_outer_EastEdge_px , FOVH_px / 2 + 1./2]

if keyword_set(printout) then begin
      print,'--------------------------------------------------------------------'
      print,'Distance [Rsun] between Sun center and:'
      print,'--------------------------------------------------------------------'
      print,'                       East-Edge       Center          West-Edge'
      print,'WISPR INNER FOV:',[D_SunCenter_inner_EastEdge,D_SunCenter_inner_Center,D_SunCenter_inner_WestEdge]/rsun
      print,'WISPR OUTER FOV:',[D_SunCenter_outer_EastEdge,D_SunCenter_outer_Center,D_SunCenter_outer_WestEdge]/rsun
      print,'--------------------------------------------------------------------'      
      print
      print,'Distance [pixels] between Sun center and:'
      print,'--------------------------------------------------------------------'
      print,'                       East-Edge       Center          West-Edge'
      print,'WISPR INNER FOV:',[D_SunCenter_inner_EastEdge_px,D_SunCenter_inner_Center_px,D_SunCenter_inner_WestEdge_px]
      print,'WISPR OUTER FOV:',[D_SunCenter_outer_EastEdge_px,D_SunCenter_outer_Center_px,D_SunCenter_outer_WestEdge_px]
      print
      print,'--------------------------------------------------------------------'
      print,'Coordinates of Sun Center in pixels in each FOV, taking (0,0) in the lower left edge of each FOV:'
      print,'Inner:',Pos_SunCenter_px_inner
      print,'Outer:',Pos_SunCenter_px_outer
endif

Distances_SUN_FOV_inner_Rsun = [D_SunCenter_inner_EastEdge,D_SunCenter_inner_Center,D_SunCenter_inner_WestEdge]/rsun
Distances_SUN_FOV_outer_Rsun = [D_SunCenter_outer_EastEdge,D_SunCenter_outer_Center,D_SunCenter_outer_WestEdge]/rsun
Distances_SUN_FOV_inner_px   = [D_SunCenter_inner_EastEdge_px,D_SunCenter_inner_Center_px,D_SunCenter_inner_WestEdge_px]
Distances_SUN_FOV_outer_px   = [D_SunCenter_outer_EastEdge_px,D_SunCenter_outer_Center_px,D_SunCenter_outer_WestEdge_px]

return
end

pro openfile,filename
openw,1,filename
return
end

pro closefile
close,/all
return
end

pro writedata,terminal=terminal,SciOrbNum=SciOrbNum,epocharray=epocharray,etarray=etarray
common output,listtype
common constants,c,rsun,au
common indexes,i,j
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun, Distances_SUN_FOV_inner_px, Distances_SUN_FOV_outer_px,px_inner_arcsec,px_outer_arcsec

if NOT keyword_set(terminal) then begin

if listtype eq 'full' then $  
  printf,1,SciOrbNum[i],epocharray[i,j],etarray[i,j],$
           sun_spp_vector_J2000, dist_SUN_SPP, dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Pos_SunCenter_px_inner, Pos_SunCenter_px_outer,$
           Distances_SUN_FOV_inner_px, Distances_SUN_FOV_outer_px,$
           Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,$
           format='(1(I4),1(A22),1(E19.10),4(F14.2),2(F10.4),2(F10.2),10(F10.1),6(F10.4))'
if listtype eq 'short' then $  
  printf,1,SciOrbNum[i],epocharray[i,j],$
           dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Distances_SUN_FOV_inner_Rsun([0,2]), Distances_SUN_FOV_outer_Rsun([0,2]),$
           format='(1(I4),1(A22),2(F10.4),2(F10.2),4(F10.4))'

endif else begin

if listtype eq 'full' then $  
  print   ,SciOrbNum,epocharray,etarray,$
           sun_spp_vector_J2000, dist_SUN_SPP, dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Pos_SunCenter_px_inner, Pos_SunCenter_px_outer,$
           Distances_SUN_FOV_inner_px, Distances_SUN_FOV_outer_px,$
           Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun,$
           format='(1(I4),1(A22),1(E19.10),4(F14.2),2(F10.4),2(F10.2),10(F10.1),6(F10.4))'
if listtype eq 'short' then $  
  print   ,SciOrbNum,epocharray,$
           dist_SUN_SPP/rsun, dist_SUN_SPP/au,$
           long_start, lat_start,$
           Distances_SUN_FOV_inner_Rsun([0,2]), Distances_SUN_FOV_outer_Rsun([0,2]),$
           format='(1(I4),1(A22),2(F10.4),2(F10.2),4(F10.4))'
endelse

return
end

pro writekeys,terminal=terminal
common output,listtype

if listtype eq 'full'  then begin  
string1=' Orb     DATE       TIME           ET            SObsJ2000x    SObsJ2000y    SObsJ2000z        D-SObs    D-SObs    D-SObs       LON       LAT    CtrI_x    CtrI_y    CtrO_x    CtrO_y      S-IE      S-IC      S-IW      S-OE      S-OC      S-OW      S-IE      S-IC      S-IW      S-OE      S-OC      S-OW'
string2='   #                                                   [km]          [km]          [km]          [km]    [Rsun]      [AU]     [deg]     [deg]      [px]      [px]      [px]      [px]      [px]      [px]      [px]      [px]      [px]      [px]    [Rsun]    [Rsun]    [Rsun]    [Rsun]    [Rsun]    [Rsun]' 
endif else begin
string1=' Orb     DATE       TIME      D-SObs    D-SObs       LON       LAT      S-IE      S-IW      S-OE      S-OW'
string2='   #                          [Rsun]      [AU]     [deg]     [deg]    [Rsun]    [Rsun]    [Rsun]    [Rsun]' 
endelse

if not keyword_set(terminal) then begin
dash_line
printf,1,string1
printf,1,string2
dash_line
endif else begin
dash_line,/terminal
print,string1
print,string2
dash_line,/terminal
endelse

return
end

pro dash_line,terminal=terminal
  common output,listtype
if listtype eq 'full'  then $
string1='-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
if listtype eq 'short' then $
string1='----------------------------------------------------------------------------------------------------------'
if not keyword_set(terminal) then printf,1,string1
if     keyword_set(terminal) then print ,  string1
return
end

pro point_line
common output,listtype
if listtype eq 'full'  then $
printf,1,'  ...........................................................................................................................................................................................................................................................................................................'
if listtype eq 'short' then $
printf,1,'  .......................................................................................................'
return
end

pro loadconstants
common constants,c,rsun,au
; Set useful constants
  c    = 299792.458      ; km/sec
  rsun = 695700.         ; km 
  au   = 149.597870700e6 ; km
return
end

FUNCTION get_dis, camera_xyz, vec
	dis = SQRT( (camera_xyz[0]-vec[0])^2 + (camera_xyz[1]-vec[1])^2 + (camera_xyz[2]-vec[2])^2 )
	RETURN,dis
END

pro Create_FITS
common constants,c,rsun,au
common spp_numbers, sun_spp_vector_J2000, dist_SUN_SPP, long_start, lat_start, Pos_SunCenter_px_inner, Pos_SunCenter_px_outer, Distances_SUN_FOV_inner_Rsun, Distances_SUN_FOV_outer_Rsun, Distances_SUN_FOV_inner_px, Distances_SUN_FOV_outer_px,px_inner_arcsec,px_outer_arcsec
common SynthFITS,hdr_Inner_0,hdr_Outer_0,img_Inner_0,img_Outer_0,datadir,epoch,et

; Create HDR equal to "TestImages" Header
hdr_Inner = hdr_Inner_0
hdr_Outer = hdr_Outer_0

; Create filenames
filename_Inner = hdr_Inner.DATE_OBS+'_Inner_Blank.fits'
filename_Outer = hdr_Inner.DATE_OBS+'_Outer_Blank.fits'

; Update .FILENAME
hdr_Inner.FILENAME = filename_Inner
hdr_Outer.FILENAME = filename_Outer

; Update .DATE_OBS
hdr_Inner.DATE_OBS =  strmid(epoch,0,10)+'T'+strmid(epoch,12,10)
hdr_Outer.DATE_OBS =  strmid(epoch,0,10)+'T'+strmid(epoch,12,10)

; Update .DSUN_OBS, distance [m] to Sun Center
hdr_Inner.DSUN_OBS = dist_SUN_SPP * 1.d3
hdr_Outer.DSUN_OBS = dist_SUN_SPP * 1.d3

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
            'CRLN_OBS',long_start,$
            'CRLT_OBS', lat_start,$
            'crota'   , 0.,$
    	    'posx_obs',sun_spp_vector_J2000[0],$
    	    'posy_obs',sun_spp_vector_J2000[1],$
    	    'posz_obs',sun_spp_vector_J2000[2]) ; km

hdr_Outer = create_struct(hdr_Outer,$
            'CRLN_OBS',long_start,$
            'CRLT_OBS', lat_start,$
            'crota'   , 0.,$
            'posx_obs',sun_spp_vector_J2000[0],$
    	    'posy_obs',sun_spp_vector_J2000[1],$
    	    'posz_obs',sun_spp_vector_J2000[2]) ; km

; Write FITS files to disk
mwritefits,hdr_Inner,img_Inner_0,outfile=datadir+filename_Inner
mwritefits,hdr_Outer,img_Outer_0,outfile=datadir+filename_Outer

stop

return
end
