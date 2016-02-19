begin

  dom_name  = "surfdata_0300x0300.nc"
  fmx_name  = "fmax_100m.nc"  

  fupdate     = True 
  ;
  in     = addfile(dom_name,"w")
  infmx  = addfile(fmx_name,"r")

  fmx      = infmx->fmax 

  in->FMAX = (/fmx/)
end
