load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;double time ( time )
;         long_name :	observation time
;         units :	days since 1972-01-01 00:00:00
;         calendar :	noleap

begin
;  in = addfile("./newforc/2009-03.nc","r")
;  tbot = in->TBOT
;  printVarSummary(tbot&time)
;  exit

 diri = "/daten01/z4/database/clm3.5/Rur_NRW/atm_forcing/"
 odir = "/daten01/z4/database/cesm/inputdata/atm/datm7/CLM1PT_data/300x300pt_NRW/"
 year = "2014"
 fils  = systemfunc("ls "+diri+year+"*")
 fregrid = False    ;WILL TAKE FOREVER USING NCL INTERPOLATIONS 
 nfils = dimsizes(fils)
 do nf = 0, nfils-1
   in = addfile(fils(nf),"r")
   t  = in->TBOT
   q  = in->QBOT
   p  = in->PSRF
   lw = in->FLDS
   sw = in->FSDS
   pr = in->PRECTmms
    w  = in->WIND
   lat = in->LONGXY
   lon = in->LATIXY  
   edgen = in->EDGEN
   edges = in->EDGES
   edgee = in->EDGEE
   edgew = in->EDGEW
   ;
   ndim = dimsizes(t)
   ntime = ndim(0)
   nlat  = ndim(1)
   nlon  = ndim(2)

   if (fregrid) then          ;fregrid = True

     lndfil   = addfile("./domain.lnd.300x300pt_NRW_navy.nc","r")
     lon_reg  = lndfil->xc(0,:)
     lat_reg  = lndfil->yc(:,0)
     lmask    = lndfil->mask
     ndim_reg = dimsizes(lmask)     

     tbot    = new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")
     wind    = new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")
     qbot    = new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")
     prectmms= new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")
     flds    = new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")
     fsds    = new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")
     psrf    = new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")
     zbot    = new((/ntime,ndim_reg(0),ndim_reg(1)/),"float")

     do it = 0, ntime-1
       print("Regridding for : " + sprinti("%0.3i",it))
       ;
       temp = rcm2rgrid(lat,lon,t(it,:,:),lat_reg,lon_reg,0)
       tbot(it,:,:) = (/temp/)
       delete(temp)
       ;
       temp = rcm2rgrid(lat,lon,w(it,:,:),lat_reg,lon_reg,0)
       wind(it,:,:) = (/temp/)
       delete(temp)
       ;
       temp = rcm2rgrid(lat,lon,q(it,:,:),lat_reg,lon_reg,0)
       qbot(it,:,:) = (/temp/)
       delete(temp)
       ;
       temp = rcm2rgrid(lat,lon,pr(it,:,:),lat_reg,lon_reg,0)
       prectmms(it,:,:) = (/temp/)
       delete(temp)
       ;
       temp = rcm2rgrid(lat,lon,sw(it,:,:),lat_reg,lon_reg,0)
       fsds(it,:,:) = (/temp/)
       delete(temp)
       ;
       temp = rcm2rgrid(lat,lon,lw(it,:,:),lat_reg,lon_reg,0)
       flds(it,:,:) = (/temp/)
       delete(temp)
       ;
       temp = rcm2rgrid(lat,lon,pr(it,:,:),lat_reg,lon_reg,0)
       psrf(it,:,:) = (/temp/)
       delete(temp)
  
       zbot(it,:,:) = 10.0
  end do

  end if      ;fregrid = True
  
 
   delete_VarAtts(lat,-1)
   delete_VarAtts(lon,-1)
   delete_VarAtts(t,-1)
   delete_VarAtts(q,-1)
   delete_VarAtts(p,-1)

   edgee@long_name = "eastern edge in atmospheric data"
   edgee@units     = "degrees_east"
   edgee@mode      = "time-invariant"
 
   edgew@long_name = "western edge in atmospheric data"
   edgew@units     = "degrees_east"
   edgew@mode      = "time-invariant"

   edgen@long_name = "northern edge in atmospheric data"
   edgen@units     = "degrees_north"
   edgen@mode      = "time-invariant"

   edges@long_name = "southern edge in atmospheric data"
   edges@units     = "degrees_north"
   edges@mode      = "time-invariant"

    lat@long_name = "latitude"
    lat@units     = "degrees_north"
    lon@long_name = "longitude"
    lon@units     = "degrees_east"

    lat@mode    = "time-invariant"
    lon@mode    = "time-invariant"
    sw@mode     = "time-dependent"
    lw@mode     = "time-dependent"
    pr@mode     = "time-dependent"
    w@mode     = "time-dependent"

    copy_VarMeta(sw, t)
    t@long_name = "temperature at the lowest atm level"
    t@units     = "K"

    copy_VarMeta(t, p)
    p@long_name = "pressure at the lowest atm level"
    p@units     = "Pa"

    copy_VarMeta(t, q)
    q@long_name = "specific humidity at the lowest atm level"
    q@units     = "kg/kg"
    
    rh = relhum(t,q,p)
    rh = rh < 100
    copy_VarMeta(t, rh)
    rh@long_name = "relative humidity at the lowest atm level"
    rh@units     = "%"

    z  = t 
    z@long_name = "height at the lowest atm level"
    z@units     = "m"
    z = 10.0

    time = new(ntime,"double")
    time@long_name = "observation time"
    time@units     = "days since " + year +"-" + sprinti("%0.2i",nf+1) + "-01 00:00:00"
    time@calendar  = "noleap"

    scal = 1.d/24.d
    if (nf.eq.0) then
      it_st = 0.d
    end if
    do it = 0, ndim(0)-1
      time(it) = it_st + it*scal
    end do
    it_st  = 0.d ;reset time(ndim(0)-1)

    print(" " + it_st)
    print(fils(nf) + " " +  time(ndim(0)-1))

    ;;
   month = nf+1
   fout = addfile(odir+year+"-"+sprinti("%0.2i",month)+".nc","c")
;
   setfileoption(fout,"DefineMode",True)
   fAtt               = True            ; assign file attributes
   fAtt@Author        = "P. Shrestha"
   fAtt@title         = "Z4, TR32 - COSMO DE_laf Regridded data for CESM"
   fAtt@source_file   = "DWD Analysis"
   fAtt@Conventions   = "CF-1.0"
   fileattdef( fout, fAtt )            ; copy file attributes
;
   dimNames = (/"scalar","time", "lat", "lon"/)
   dimSizes = (/ 1, ntime  ,  nlat,  nlon /)
   dimUnlim = (/ False, True , False, False/)
   filedimdef(fout,dimNames,dimSizes,dimUnlim)
;
   filevardef(fout, "time", typeof(time), "time")
   filevardef(fout, "EDGEE", typeof(edgee), "scalar")
   filevardef(fout, "EDGEW", typeof(edgew), "scalar")
   filevardef(fout, "EDGEN", typeof(edgen), "scalar")
   filevardef(fout, "EDGES", typeof(edges), "scalar")

   filevardef(fout, "LONGXY" ,typeof(lon),(/"lat","lon"/))
   filevardef(fout, "LATIXY" ,typeof(lat),(/"lat","lon"/))

   filevardef(fout, "TBOT" ,typeof(t),(/"time","lat","lon"/))
   filevardef(fout, "WIND" ,typeof(w),(/"time","lat","lon"/))
   filevardef(fout, "RH" ,typeof(rh),(/"time","lat","lon"/))
   filevardef(fout, "PRECTmms" ,typeof(pr),(/"time","lat","lon"/))
   filevardef(fout, "FSDS" ,typeof(sw),(/"time","lat","lon"/))
   filevardef(fout, "FLDS" ,typeof(lw),(/"time","lat","lon"/))
   filevardef(fout, "PSRF" ,typeof(p),(/"time","lat","lon"/))
   filevardef(fout, "ZBOT" ,typeof(z),(/"time","lat","lon"/))
 ;
   filevarattdef(fout, "LONGXY"  , lon)
   filevarattdef(fout, "LATIXY"  , lat)
   filevarattdef(fout, "TBOT"  , t)
   filevarattdef(fout, "WIND"  , w)
   filevarattdef(fout, "RH"  , rh)
   filevarattdef(fout, "PRECTmms"  , pr)
   filevarattdef(fout, "FSDS"  , sw)
   filevarattdef(fout, "FLDS"  , lw)
   filevarattdef(fout, "PSRF"  , p)
   filevarattdef(fout, "ZBOT"  , z)
   filevarattdef(fout, "time"  , time)
   filevarattdef(fout, "EDGEE"  , edgee)
   filevarattdef(fout, "EDGEW"  , edgew)
   filevarattdef(fout, "EDGEN"  , edgen)
   filevarattdef(fout, "EDGES"  , edges)

   setfileoption(fout,"DefineMode",False)


   fout->EDGEN    =(/edgen/)
   fout->EDGEE    =(/edgee/)
   fout->EDGES    =(/edges/)
   fout->EDGEW    =(/edgew/)
   fout->time     = (/time/)
   fout->LONGXY   =(/lon/)
   fout->LATIXY   =(/lat/)
   fout->TBOT     =(/t/)
   fout->RH       = (/rh/)
   fout->WIND     =(/w/)
   fout->PRECTmms =(/pr/)
   fout->FSDS     =(/sw/)
   fout->FLDS     =(/lw/)
   fout->PSRF     =(/p/)
   fout->ZBOT     = (/z/)
;;
   delete([/z,in, t,p,q,rh,sw,lw,pr,w,lon,lat, ndim, time/])
  end do
end
