load "Func/plot_map_res.ncl"

begin
def = "BS"
   figtype = "oldps"
   figname = "Pic/"+def+"/kinetic_energy_conver"

  ff = addfile("/work/L.xhm/anal_mean/Comp/Data/"+def+"_kinetic_energy_conver_DJFM.nc", "r")
u2_dU_bardx   = dim_avg_n_Wrap(ff->u2_dU_bardx ,0) *10^4
uv_dU_bardy   = dim_avg_n_Wrap(ff->uv_dU_bardy ,0) *10^4
uv_dV_bardx   = dim_avg_n_Wrap(ff->uv_dV_bardx ,0) *10^4
v2_dV_bardy   = dim_avg_n_Wrap(ff->v2_dV_bardy ,0) *10^4
kp            = dim_avg_n_Wrap(ff->kp ,0)

copy_VarCoords(kp, u2_dU_bardx)
copy_VarCoords(kp, uv_dU_bardy)
copy_VarCoords(kp, uv_dV_bardx)
copy_VarCoords(kp, v2_dV_bardy)
 
       kp  = kp*10^4
copy_VarCoords(u2_dU_bardx, kp)

; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  ;mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = "Dec" ;(/"JJA", "SON", "DJF", "MAM"/)
  plot = new((/6,dimsizes(mm)/), graphic)

  Levels             = (/-6., -5., -4., -3., -2., -1.5, -1., -.5, \
                          .5,  1., 1.5,  2.,  3.,   4.,  5.,  6./)
  res  = PLOT_map_res(Levels, 8, -10., 75., 80., 240.) 
  ;res  = PLOT_map_res(Levels, 8, -10., 75., 30., 300.) 
  res2 = res
  res2@lbLabelBarOn         = False


    do i = 0,0 ;dimsizes(mm)-1
;       ii = i+7
      res2@gsnLeftString    = "(a) RHS term 1" 
       plot(0,i)  = gsn_csm_contour_map_ce(wks, u2_dU_bardx(:,:), res2)

      res2@gsnLeftString    = "(b) RHS term 2" 
       plot(1,i)  = gsn_csm_contour_map_ce(wks, uv_dU_bardy(:,:), res2)

      res2@gsnLeftString    = "(c) RHS term 3" 
       plot(2,i)  = gsn_csm_contour_map_ce(wks, uv_dV_bardx(:,:), res2)

      res2@gsnLeftString    = "(d) RHS term 4" 
       plot(3,i)  = gsn_csm_contour_map_ce(wks, v2_dV_bardy(:,:), res2)

      res@gsnLeftString     = "(a) 200hPa Barotropic Conversion (m~S~2~N~/s~S~3~N~)"
       plot(5,i)  = gsn_csm_contour_map_ce(wks, kp(:,:), res)
    end do

  resP = True
  gsn_panel(wks, plot(5,0), (/1,1/), resP)

  resP@gsnPanelLabelBar     = True
  resP@lbLabelFontHeightF   = 0.01               ; make labels smaller
  gsn_panel(wks, plot(:3,0), (/2,2/), resP)


;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end
