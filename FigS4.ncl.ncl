load "Func/plot_map_res.ncl"

begin
def = "BS"
   figtype = "oldps"
   figname = "Pic/"+def+"/wind_shear"

  ff = addfile("/work/L.xhm/anal_mean/Comp/Data/"+def+"_kinetic_energy_conver_DJFM.nc", "r")
dU_bardx   = dim_avg_n_Wrap(ff->dU_bardx ,0) *10^6
dU_bardy   = dim_avg_n_Wrap(ff->dU_bardy ,0) *10^6
dV_bardx   = dim_avg_n_Wrap(ff->dV_bardx ,0) *10^6
dV_bardy   = dim_avg_n_Wrap(ff->dV_bardy ,0)

copy_VarCoords(dV_bardy, dU_bardx)
copy_VarCoords(dV_bardy, dU_bardy)
copy_VarCoords(dV_bardy, dV_bardx)

dV_bardy   = dV_bardy *10^6
copy_VarCoords(dV_bardx, dV_bardy)

 

; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  ;mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = "Dec" ;(/"JJA", "SON", "DJF", "MAM"/)
  plot = new((/6,dimsizes(mm)/), graphic)

  Levels             = (/-6., -5., -4., -3., -2., -1.5, -1., -.5, \
                          .5,  1., 1.5,  2.,  3.,   4.,  5.,  6./)*2.
  res  = PLOT_map_res(Levels, 8, -10., 75., 80., 240.) 
  ;res  = PLOT_map_res(Levels, 8, -10., 75., 30., 300.) 
  res2 = res
  res2@lbLabelBarOn         = False


    do i = 0,0 ;dimsizes(mm)-1
;       ii = i+7
      res2@gsnLeftString    = "(a) dUdx" 
       plot(0,i)  = gsn_csm_contour_map_ce(wks, dU_bardx(:,:), res2)

      res2@gsnLeftString    = "(b) dUdy"
       plot(1,i)  = gsn_csm_contour_map_ce(wks, dU_bardy(:,:), res2)

      res2@gsnLeftString    = "(c) dVdx"
       plot(2,i)  = gsn_csm_contour_map_ce(wks, dV_bardx(:,:), res2)

      res2@gsnLeftString    = "(d) dVdy"
       plot(3,i)  = gsn_csm_contour_map_ce(wks, dV_bardy(:,:), res2)

    end do

  resP = True

  resP@gsnPanelLabelBar     = True
  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller
  gsn_panel(wks, plot(:3,0), (/2,2/), resP)


;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end
