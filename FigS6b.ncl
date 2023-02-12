; Local Walker Circulation Ref:
; https://doi.org/10.1175/JCLI-D-16-0525.1
load "plot_map_res.ncl"

undef ("CompData")
function CompData(case)
local data
begin
   pai  = 3.1415926
   pwd  = "/work/L.xhm/anal_mean/Comp/Data/"
   fi   = addfile(pwd+case+"div_u_DJFM.nc", "r")
   lev  = fi->lev
   lat  = fi->lat
   lon  = fi->lon  
   data = fi->udiv(:, :, :, {-10:10}, :)
   res_data = dim_avg_n_Wrap(data,3)

   sum_data = new(dimsizes(res_data), typeof(res_data))
   dp   = dpres_plevel(lev, max(lev), min(lev), 0)*100.  ; hPa->Pa 

   do k = 0,dimsizes(lev)-1
      do y = 0,dimsizes(lon)-1
         para = tofloat(2. *pai *6371.*10.^3  /9.8)
         sum_data(:,:,k,y) = dim_sum_wgt_n(res_data(:,:,k:,y), dp(k:),  1, 2)*para
      end do
   end do
   
   sum_data = sum_data*10.^-10
   sum_data!2 = "lev"
   sum_data!3 = "lon"
   sum_data&lev = lev
   sum_data&lon = lon
   return( sum_data )
end

undef("SigTest")
function SigTest(var_ctl, var_def)
local var_ctl,var_def,dim,dim_name 
begin
   var_prob = ttest(dim_avg_n(var_def,0), dim_variance_n(var_def,0), dimsizes(var_def(:,0,0,0)),\
                    dim_avg_n(var_ctl,0), dim_variance_n(var_ctl,0), dimsizes(var_ctl(:,0,0,0)), True, False)
   var_prob@_FillValue = default_fillvalue("float")

   var_prob1  = var_prob
   var_prob2  = var_prob
   var_prob1  = where(var_prob.lt.0.05, var_prob, var_prob@_FillValue)
   var_prob2  = where(var_prob.lt.0.1,  var_prob, var_prob@_FillValue)

   var_dif     = dim_avg_n(var_def,0) - dim_avg_n(var_ctl,0)
   var_ctl_avg = dim_avg_n(var_ctl,0)
      
   rt_data = (/var_ctl_avg, var_dif, var_prob1, var_prob2/)
   copy_VarCoords(var_ctl(0,0,:,:), rt_data(0,0,:,:))
   return(rt_data)
end

   
undef("PLOT")
procedure PLOT(ctl, def, name, yyB, yyE)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = name

   ff_name = name+".nc"
   if(.not.fileexists(ff_name))then
      ;; Stream Func
      data_ctl = CompData(ctl)
      data_def = CompData(def)
      data1 = SigTest(data_ctl, data_def)

      ff_op = addfile(ff_name, "c")
      ff_op->data1 = data1
   else
      ff_in = addfile(ff_name, "r")
      data1 = ff_in->data1
   end if


plot_data1 = data1(1,:,:,:)
plot_data1 = where(data1(3,:,:,:).lt.0.1,   plot_data1, data1@_FillValue)

lev = data1&lev
lon = data1&lon

; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  ;mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = "DJFM"
  plot  = new(dimsizes(mm), graphic)
  plotT = new(dimsizes(data1(1,:,:,:)), graphic)


    res                      = True             ; plot mods desired
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnAddCyclic         = False
    res@gsnLeftString        = ""
    res@tiMainFontHeightF    = 0.02
    res@tiYAxisString        = "Depth(m)"

    res@cnFillOn             = True               ; turn on color for contours
;    res@cnFillMode           = "CellFill"
    res@cnLinesOn            = False              ; turn off contour lines
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    res@cnInfoLabelOn        = False

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnFillPalette        = ReadColorMap2("NCV_jet", 8)

;    res@lbLabelBarOn         = False
    res@lbLabelFontHeightF   = 0.01              ; make labels smaller
    res@vpWidthF             = 0.58
    res@vpHeightF            = 0.4

    res@tmXTOn               = False
    res@tmYROn               = False
    res@tmBorderThicknessF   = .8
    res@tmXBMajorThicknessF  = 1.2
    res@tmYLMajorThicknessF  = 1.2
    res@tmXBLabelFontHeightF = 0.012
    res@tmYLLabelFontHeightF = 0.012
                                                   
    res@trXMinF              = 30.
    res@trXMaxF              = 270.
    res@trYMinF              = 100.
    res@trYReverse           = True
;    res@gsnRightString       = "Pa/s"
    res@tiYAxisString        = "Height(hPa)"

    res2              = True                         ; plot mods desired
    res2@gsnDraw      = False                        ; don't draw yet
    res2@gsnFrame     = False                        ; don't advance frame yet
    res2@gsnLeftString= ""
    res2@gsnContourZeroLineThicknessF = 2.     ; doubles thickness of zero contour
    res2@gsnContourNegLineDashPattern = 11
    res2@gsnContourPosLineDashPattern = 0
    res2@cnLevelSelectionMode = "ExplicitLevels"
    res2@cnLevels             = ispan(-15,15,3)
    res2@cnLineLabelBackgroundColor = -1
    res2@cnLineColor          = "black"
    res2@cnInfoLabelOn        = False
    res2@cnInfoLabelFontColor = "black"
    res2@cnLineLabelFontColor = "black"

  RES3 = PLOT_sig_res("grey66", .009)
  resdot = RES3
  RES3 = PLOT_sig_res("grey66", .01)
  resdot2= RES3


  resP = True
;  resP@gsnPanelLabelBar     = True
;  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller

    do i = 0,dimsizes(mm)-1
       ii = i
       res@cnLevels  = (/-2., -1.6, -1.2, -1., -.7, -.5, -.3, -.1, \
                          .1,   .3,   .5,  .7,  1., 1.2, 1.6,  2./)
                       

      res@gsnLeftString     = "(b) Local Walker Circulation"
       plot(i)  = gsn_csm_contour(wks, plot_data1(ii,:,:), res)
       cli      = gsn_csm_contour(wks,   data1(0,ii,:,:), res2)
       overlay(plot(i), cli)
       delete([/cli/])

    end do

    gsn_panel(wks, (/plot(0)/),  (/1,1/),  resP)

;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end


begin

   ctl = "B2000_CAM5_CTR_"
   yyB = "0011"
   yyE = "0030"

   def = "B2000_CAM5_MC2BS_"
PLOT(ctl, def, "FigS6b", yyB, yyE)


end
