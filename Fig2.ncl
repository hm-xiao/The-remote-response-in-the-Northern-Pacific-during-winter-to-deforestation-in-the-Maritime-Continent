load "plot_map_res.ncl" 
load "kinetic_energy_conv_func.ncl" 

undef("ReadData")
function ReadData(area, name, lev_cc, var)
begin
   fin = addfile("/work/home/L.xhm/LBM_standard/Out/"+area+"/"+name+".nc", "r")
   lev = fin->lev
   levID = ind(lev.eq.lev_cc)
   if(var.eq."p")then
      data  = rm_single_dims(fin->$var$)
   else
      data  = rm_single_dims(fin->$var$(:,levID,:,:))
      aa    = data
      if(var.eq."psi")then
         aa = data*10^-6
      else if(var.eq."chi")then
         aa = data*10^-5
      end if
      end if
      copy_VarCoords(data,aa)
      data = aa
      delete(aa)
   end if
   return(data)
end


undef("PLOT")
procedure PLOT(area, bs_file, lev_cc, FigName, Levels_contour, Levels_shading)
begin
 name = area+"_t42l20_"+bs_file
ffile = "/work/home/L.xhm/LBM_standard/Out/"+area+"/"+name+".nc" 

   ff_name = FigName+".nc"
   if(.not.fileexists(ff_name))then

      fin = addfile(ffile, "r")
      lev = fin->lev
      levID = ind(lev.eq.lev_cc)

      uu = fin->u(:,levID,::-1,:)
      vv = fin->v(:,levID,::-1,:)

      ff2 = addfile("/work/home/L.xhm/LBM_standard/Bs/CESM."+bs_file+".t42.nc","r")
      u_ctl = dim_avg_n_Wrap(ff2->u(6:9,levID,:,:),0)
      v_ctl = dim_avg_n_Wrap(ff2->v(6:9,levID,:,:),0)

      ff_op = addfile(ff_name, "c")
      ff_op->uu=uu
      ff_op->vv=vv
      ff_op->u_ctl=u_ctl
      ff_op->v_ctl=v_ctl
   else
      ff_in = addfile(ff_name, "r")
      uu = ff_in->uu
      vv = ff_in->vv
      u_ctl = ff_in->u_ctl
      v_ctl = ff_in->v_ctl
   end if

   lat = tofloat(u_ctl&lat)
   lon = tofloat(u_ctl&lon)

;=============================
; POT
;=============================
   figtype  = "oldps" 
   figname  = FigName
    wks  = gsn_open_wks(figtype,figname)

   resct = PLOT_contour_res(Levels_contour,  "black")
   res   = PLOT_map_res(Levels_shading, 8,  -10., 75., 80., 240.)

   res@lbLabelBarOn       = False            ; turn off individual cb's

   delete(res@cnFillPalette)
   ;res@cnFillPalette        = ReadColorMap2("NCV_jet", 9)
   res@cnFillPalette        = ReadColorMap2("BlWhRe", 9)

  data_kpe = new(dimsizes(uu), typeof(uu))
    do i = 0,dimsizes(uu(:,0,0))-1
       data_kpe(i,:,:) = KPE(u_ctl, v_ctl, uu(i,:,:), vv(i,:,:)) 
    end do
   
    data_kpe = data_kpe *10^5
    copy_VarCoords(uu, data_kpe)

    sf  =  ReadData(area, name, lev_cc, "psi")
    sf  = dim_rmvmean_n_Wrap(sf,2)

;;;;  Average PLOT------------------------------
    res@gsnLeftString        = " Avg: Day 20-49"
    plot = gsn_csm_contour_map_ce(wks, dim_avg_n_Wrap(sf(29:58,:,:),0), res)
     aa := gsn_csm_contour(wks, dim_avg_n_Wrap(data_kpe(29:58,:,:),0), resct)
     overlay(plot, aa)

;    draw(plot)
;    frame(wks)


;;;;  Panel PLOT------------------------------
   plot2 = new(9, graphic)
char = (/"a","b","c","d","e","f","g","h","i"/)

   j = 0
   do i = 1,9,2
      res@gsnLeftString   = "("+char(j)+") Day "+(i+1)
      plot2(j) = gsn_csm_contour_map_ce(wks, sf(i,:,:), res)
      aa := gsn_csm_contour(wks, data_kpe(i,:,:), resct)
      overlay(plot2(j), aa)
      j = j+1
   end do

   do i = 14,24,5
      res@gsnLeftString   = "("+char(j)+") Day "+(i+1)
      plot2(j) = gsn_csm_contour_map_ce(wks, sf(i,:,:), res)
      aa := gsn_csm_contour(wks, data_kpe(i,:,:), resct)
      overlay(plot2(j), aa)
      j = j+1
   end do


   resP                    = True                 ; modify the panel plot
   resP@gsnPanelLabelBar   = True                 ; add common colorbar
   gsn_panel(wks, (/plot2(0),plot2(1),plot2(2),plot2(3),plot2(4),plot2(5),\
                    plot2(6),plot2(7),plot/),(/3,3/),resP)               ; now draw as one plot


;========================Transform Image File========================
 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end

begin
;procedure PLOT(area, bs_file, lev_cc, FigName, Levels_contour, Levels_shading)
  Levels_contour = (/-4., -2., -1., -.6, -.4, -.2, .2, .4, .6, 1.,  2.,  4./)
  Levels_shading = (/-6., -5., -4., -3., -2., -1.5, -1., -.5, \
                      .5,  1., 1.5,  2.,  3.,   4.,  5.,  6./)*.2


fr = (/"MC"/)
bs = (/"cli"/)

  PLOT(fr,  bs, 200, "Fig2", Levels_contour*2.,    Levels_shading)
end
     
