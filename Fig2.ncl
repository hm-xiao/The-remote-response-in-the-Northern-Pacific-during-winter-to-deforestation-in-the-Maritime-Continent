load "/work/L.xhm/anal_ensemble/Func/plot_map_res.ncl"

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
procedure PLOT(area, bs_file, lev_cc, LeftString, var1,var2,var3,var4, Levels, RefMagnitudeF, FigName)
begin
   name = area+bs_file+"_t42l20_60dd"
   system("mkdir -p Pic/"+name) ;+"/png/")

   data1  = ReadData(area, name, lev_cc, var1)
   data2  = ReadData(area, name, lev_cc, var2)
   data3  = ReadData(area, name, lev_cc, var3)
   data4  = ReadData(area, name, lev_cc, var4)

   data1  = dim_rmvmean_n_Wrap(data1,2)

;=============================
; POT
;=============================
   figtype = "oldps"; "eps"
   figname = "Pic/"+name+"/"+FigName
    wks = gsn_open_wks(figtype,figname)

;  res = PLOT_map_res(Levels, 8,  -65., 65., 0., 360.)
  res = PLOT_map_res(Levels, 8,  -10., 75., 80., 240.)
  res@lbLabelBarOn       = False            ; turn off individual cb's

;  resvc = PLOT_vc_res(2., 0.015, 0.04, RefMagnitudeF, .08)


;;;;  Average PLOT------------------------------
    res@gsnLeftString        =  "(i) Day 30-59 Avg"
       plot = gsn_csm_contour_map_ce(wks, dim_avg_n_Wrap(data1(29:58,:,:),0), res)
;       aa = gsn_csm_contour(wks, dim_avg_n_Wrap(data4(29:58,:,:),0), resct)
;       bb = gsn_csm_vector( wks, dim_avg_n_Wrap(data2(29:58,:,:),0),  dim_avg_n_Wrap(data3(29:58,:,:),0), resvc)

;       overlay(plot, aa)
;if(var1.eq."temp")then
;       overlay(plot, bb)
;end if
;      draw(plot)
;      frame(wks)


;;;;  Panel PLOT------------------------------
char = (/"a","b","c","d","e","f","g","h","i"/)
   plot2 = new(9, graphic)
   j = 0
   do i = 1,9,2
      res@gsnLeftString   = "("+char(j)+") Day "+(i+1)
      plot2(j) = gsn_csm_contour_map_ce(wks, data1(i,:,:), res)
      j = j+1
   end do

   do i = 14,24,5
      res@gsnLeftString   = "("+char(j)+") Day "+(i+1)
      plot2(j) = gsn_csm_contour_map_ce(wks, data1(i,:,:), res)
      j = j+1
   end do

   resP                    = True                 ; modify the panel plot
   resP@gsnPanelLabelBar   = True                 ; add common colorbar
   gsn_panel(wks, (/plot2(0),plot2(1),plot2(2),plot2(3),\
                    plot2(4),plot2(5),plot2(6),plot2(7),plot/),(/3,3/),resP)               ; now draw as one plot


;;;;  Solo PLOT------------------------------
;   figtype2 = "png" 
;   figname2 = "Pic/"+name+"/png/"+FigName
;    wks2 = gsn_open_wks(figtype2,figname2)

;   do i = 0,50
;      res@gsnLeftString   = "day"+(i+1)
;      plot = gsn_csm_contour_map_ce(wks2, data1(i,:,:), res)
;      aa = gsn_csm_contour(wks2, data4(i,:,:), resct)
;      bb = gsn_csm_vector( wks2, data2(i,:,:),  data3(i,:,:), resvc)
;      overlay(plot, aa)
;if(var1.eq."temp")then
;     overlay(plot, bb)
;end if
;      draw(plot)
;      frame(wks2)
;   end do


;========================Transform Image File========================
 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if
end

begin
;procedure PLOT(area, bs_file, lev_cc, var1,var2,var3, Levels, FigName)
event = (/"cli"/) ;, "nina", "nino"/)
area = (/"MC_CESM_ocn_Dec"/) ;(/"MC",  "MCL", "MCO", "CP", "EP"/)
mm = (/"Dec"/) ;(/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb","Mar", "Apr", "May"/)

Levels_slp =  (/-.5, -.4, -.3, -.25, -.2, -.15, -.1, -.05,\
                .05,  .1, .15,   .2, .25,   .3,  .4,   .5/)

Levels_temp =  (/-3.6, -2.8, -2., -1.6, -1.2, -.8, -.4, -.2,\
                   .2,   .4,  .8,  1.2,  1.6,  2., 2.8, 3.6/)/4.

Levels_200z = (/ -3.,-2.,-1.5, -1.,-.8, -.6, -.4, -.2,\
                  .2, .4,  .6,  .8, 1., 1.5,  2.,   3./) *4.

 
Levels_200sf = (/ -1., -.8, -.6, -.5, -.4, -.3, -.2, -.1,\
                   .1,  .2,  .3,  .4,  .5,  .6,  .8,  1./) 
 

do i = 0,dimsizes(event)-1
;do j = 0,dimsizes(area)-1
do k = 0,dimsizes(mm)-1
   PLOT("MC_CESM_lnd_DJFM",   "", 200, "", "psi","u","v","z", Levels_200sf, 5., "200sf" )
   PLOT("MC2C4_CESM_lnd_DJFM",   "", 200, "", "psi","u","v","z", Levels_200sf/2., 5., "200sf" )
;   PLOT("MC_CESM_lnd_Dec",   "", 200, "",   "z","u","v","z", Levels_200z, 5.,  "200z" )
;   PLOT("MC_CESM_lnd_Dec",   "", 850, "", "psi","u","v","z", Levels_200sf, 5., "850sf" )
;   PLOT("MC_CESM_lnd_Dec",   "", 850, "",   "z","u","v","z", Levels_200z/4., 5.,  "850z" )
end do
;end do
end do


end
     
