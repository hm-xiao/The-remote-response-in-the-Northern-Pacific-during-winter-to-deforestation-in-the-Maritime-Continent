load "plot_map_res.ncl"

begin
   figtype = "oldps"
   figname = "FigS4"

   ff_name = figname+".nc"
   if(.not.fileexists(ff_name))then


      ff = addfile("/work/L.xhm/anal_mean/Comp/Data/BS_kinetic_energy_conver_DJFM.nc", "r")
      dU_bardx   = dim_avg_n_Wrap(ff->dU_bardx ,0) *10^6
      dU_bardy   = dim_avg_n_Wrap(ff->dU_bardy ,0) *10^6
      dV_bardx   = dim_avg_n_Wrap(ff->dV_bardx ,0) *10^6
      dV_bardy   = dim_avg_n_Wrap(ff->dV_bardy ,0)

      copy_VarCoords(dV_bardy, dU_bardx)
      copy_VarCoords(dV_bardy, dU_bardy)
      copy_VarCoords(dV_bardy, dV_bardx)

      dV_bardy   = dV_bardy *10^6
      copy_VarCoords(dV_bardx, dV_bardy)

      ff_op = addfile(ff_name, "c")
      ff_op->dU_bardx=dU_bardx
      ff_op->dU_bardy=dU_bardy
      ff_op->dV_bardx=dV_bardx
      ff_op->dV_bardy=dV_bardy
   else
      ff_in = addfile(ff_name, "r")
      dU_bardx = ff_in->dU_bardx
      dU_bardy = ff_in->dU_bardy
      dV_bardx = ff_in->dV_bardx
      dV_bardy = ff_in->dV_bardy
   end if


; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  plot = new(4, graphic)

  Levels             = (/-6., -5., -4., -3., -2., -1.5, -1., -.5, \
                          .5,  1., 1.5,  2.,  3.,   4.,  5.,  6./)*2.
  res  = PLOT_map_res(Levels, 8, -10., 75., 80., 240.) 
  res2 = res
  res2@lbLabelBarOn         = False

      res2@gsnLeftString    = "(a) dUdx" 
       plot(0)  = gsn_csm_contour_map_ce(wks, dU_bardx, res2)

      res2@gsnLeftString    = "(b) dUdy"
       plot(1)  = gsn_csm_contour_map_ce(wks, dU_bardy, res2)

      res2@gsnLeftString    = "(c) dVdx"
       plot(2)  = gsn_csm_contour_map_ce(wks, dV_bardx, res2)

      res2@gsnLeftString    = "(d) dVdy"
       plot(3)  = gsn_csm_contour_map_ce(wks, dV_bardy, res2)


  resP = True

  resP@gsnPanelLabelBar     = True
  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller
  gsn_panel(wks, plot, (/2,2/), resP)


;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end
