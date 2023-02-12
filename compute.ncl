undef("AreaMean")
function AreaMean(sst)
local rad,lat,wgtlat
begin
  rad  = 4.*atan(1.)/180.
  lat  = sst&lat
  wgtlat   = cos(lat*rad)
  sst_id   = wgt_areaave(sst, wgtlat, 1.0 ,0)
  return(sst_id)
end
