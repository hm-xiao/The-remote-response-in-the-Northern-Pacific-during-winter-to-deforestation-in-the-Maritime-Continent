undef("AboveTerrain")
function AboveTerrain(case, yyB, yyE, data_in, lev_in)
local fis, fin, ps, dim_name
begin
   pwd = "/work/L.xhm/anal_mean/Comp/Data/"
   yyBI = toint(yyB)
   yyEI = toint(yyE)

   fis = systemfunc("ls "+pwd+case+"*/cam_0001_00?0.nc")
   fin = addfiles(fis, "r")
   ListSetType (fin, "join")

   ps = fin[:]->PS(:, yyBI*12-12:yyEI*12-1, : ,:)
   ps_hPa = ps/100.
   copy_VarCoords(ps, ps_hPa)
   ps = ps_hPa

   ret_data = data_in
   res_lev  = data_in

   if(lev_in .ne. "sfc")then
   res_lev  = tofloat(lev_in)
   if(dimsizes(lev_in).gt.1)then
      do i = 0,dimsizes(lev_in)-1
         ret_data(:,:,i,:,:) = where(res_lev(i,:,:).lt.ps, ret_data(:,:,i,:,:), ret_data@_FillValue)
      end do
   else
      ret_data = where(res_lev.lt.ps, ret_data, ret_data@_FillValue)
   end if
   copy_VarCoords(data_in, ret_data)
   end if
   
   return(ret_data)
end
