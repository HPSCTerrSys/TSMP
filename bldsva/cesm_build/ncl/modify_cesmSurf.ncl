begin

  dom_name  = "surfdata_0300x0300.nc"
  fmx_name  = "fmax_100m.nc"  
  base_name = "surfdata_300x300pt_NRW_simyr2000.nc"
  fout_name = "new.nc"

  fupdate     = True 
  ;
  in     = addfile(dom_name,"r")
  inBASE = addfile(base_name,"r")
  infmx  = addfile(fmx_name,"r")

  ;names     = getfilevarnames(inBASE) 
 
; First Create a base file of your dimension
  print("rm "+ fout_name)
  system("rm "+ fout_name)
  ;print("ncks -d lsmlat,59,358 -d lsmlon,0,299 " +base_name + " -O " + fout_name)
  ;system("ncks -d lsmlat,59,358 -d lsmlon,0,299 " +base_name + " -O " + fout_name)
  ; 
  system("cp " + base_name + " " + fout_name)

; Extract existing surface data
  l_pft    = todouble(in->LANDFRAC_PFT) 
  p_msk    = tointeger(in->PFTDATA_MASK)
  s_col    = tointeger(in->SOIL_COLOR)
  p_snd    = todouble(in->PCT_SAND)
  p_cly    = todouble(in->PCT_CLAY)
  p_wtl    = todouble(in->PCT_WETLAND)
  p_lak    = todouble(in->PCT_LAKE)
  p_gla    = todouble(in->PCT_GLACIER) 
  fmx      = todouble(infmx->fmax) 
  p_urb    = todouble(in->PCT_URBAN) 
  p_pft    = todouble(in->PCT_PFT) 
  m_lai    = todouble(in->MONTHLY_LAI) 
  m_sai    = todouble(in->MONTHLY_SAI) 
  m_ht     = todouble(in->MONTHLY_HEIGHT_TOP) 
  m_hb     = todouble(in->MONTHLY_HEIGHT_BOT) 

  ;Change bare soil to urban area
  p_urb   = p_pft(16,:,:)
  p_pft(16,:,:) = 0.d

  if (fupdate) then
  out    = addfile(fout_name,"w")
  out->LANDFRAC_PFT =(/l_pft/)
  out->PFTDATA_MASK =(/p_msk/)
  out->SOIL_COLOR =(/s_col/)
  out->PCT_SAND =(/p_snd/)
  out->PCT_CLAY =(/p_cly/)
  out->PCT_WETLAND =(/p_wtl/)
  out->PCT_LAKE =(/p_lak/)
  out->PCT_GLACIER =(/p_gla/)
  out->FMAX =(/fmx/)
  out->PCT_URBAN =(/p_urb/)
  out->PCT_PFT =(/p_pft/)
  out->MONTHLY_LAI =(/m_lai/)
  out->MONTHLY_SAI =(/m_sai/)
  out->MONTHLY_HEIGHT_TOP =(/m_ht/)
  out->MONTHLY_HEIGHT_BOT =(/m_hb/)

  end if
end
