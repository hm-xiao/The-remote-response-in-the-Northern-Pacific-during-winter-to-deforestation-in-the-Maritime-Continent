load "plot_map_res.ncl"

begin
   figtype = "oldps"
   figname = "Fig3a_S3" 

   ff_name = figname+".nc"
   if(.not.fileexists(ff_name))then

      ff = addfile("/work/L.xhm/anal_mean/Comp/Data/BS_kinetic_energy_conver_DJFM.nc", "r")
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

      ff_op = addfile(ff_name, "c")
      ff_op->u2_dU_bardx=u2_dU_bardx
      ff_op->uv_dU_bardy=uv_dU_bardy
      ff_op->uv_dV_bardx=uv_dV_bardx
      ff_op->v2_dV_bardy=v2_dV_bardy
      ff_op->kp=kp
   else
      ff_in = addfile(ff_name, "r")
      u2_dU_bardx = ff_in->u2_dU_bardx
      uv_dU_bardy = ff_in->uv_dU_bardy
      uv_dV_bardx = ff_in->uv_dV_bardx
      v2_dV_bardy = ff_in->v2_dV_bardy
      kp = ff_in->kp
   end if


; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  plot = new(6, graphic)

  Levels             = (/-6., -5., -4., -3., -2., -1.5, -1., -.5, \
                          .5,  1., 1.5,  2.,  3.,   4.,  5.,  6./)
  res  = PLOT_map_res(Levels, 8, -10., 75., 80., 240.) 
  res2 = res
  res2@lbLabelBarOn         = False


;       ii = i+7
      res2@gsnLeftString    = "(a) RHS term 1" 
       plot(0)  = gsn_csm_contour_map_ce(wks, u2_dU_bardx, res2)

      res2@gsnLeftString    = "(b) RHS term 2" 
       plot(1)  = gsn_csm_contour_map_ce(wks, uv_dU_bardy, res2)

      res2@gsnLeftString    = "(c) RHS term 3" 
       plot(2)  = gsn_csm_contour_map_ce(wks, uv_dV_bardx, res2)

      res2@gsnLeftString    = "(d) RHS term 4" 
       plot(3)  = gsn_csm_contour_map_ce(wks, v2_dV_bardy, res2)

      res@gsnLeftString     = "(a) 200hPa Barotropic Conversion (m~S~2~N~/s~S~3~N~)"
       plot(4)  = gsn_csm_contour_map_ce(wks, kp, res)

  resP = True
  gsn_panel(wks, plot(4), (/1,1/), resP)

  resP@gsnPanelLabelBar     = True
  resP@lbLabelFontHeightF   = 0.01               ; make labels smaller
  gsn_panel(wks, plot(:3), (/2,2/), resP)


;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end
