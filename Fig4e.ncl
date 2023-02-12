load "plot_map_res.ncl"
load "compute.ncl"

undef ("ReadData")
function ReadData(case, var, yyB, yyE)
local Data,data,data_dt,data_avg,dim,dim_name
begin
   pwd = "/work/L.xhm/anal_mean/Comp/Data/"
   yyBI = toint(yyB)
   yyEI = toint(yyE)

   fis = systemfunc("ls "+pwd+case+"*/"+var+"_0001_00?0.nc")
   fin = addfiles(fis, "r")
   ListSetType (fin, "join")

   Data = fin[:]->$var$(:, yyBI*12-12:yyEI*12-1, :, {50:75}, {150:200})
   dim  = dimsizes(Data)
   dim_name = getvardims(Data)

   Data@_FillValue     = default_fillvalue("float")
   Data@missing_value  = default_fillvalue("float")

   data = reshape(Data, (/dim(0)*dim(1), dim(2), dim(3), dim(4)/))
   data = runave_n(data, 4, 0, 0)

   res_data  = reshape(data(6-1:dim(0)*dim(1)-7-1,:,:,:), (/dim(0)*dim(1)/12-1, 12, dim(2),dim(3),dim(4)/)) ;  #####
   res_data!2 = dim_name(2)
   res_data!3 = dim_name(3)
   res_data!4 = dim_name(4)
   res_data&$dim_name(2)$ = Data&$dim_name(2)$
   res_data&$dim_name(3)$ = Data&$dim_name(3)$
   res_data&$dim_name(4)$ = Data&$dim_name(4)$
   delete([/Data, data/])
   Data = res_data(:,7:7,:,:,:)

   fi_lnd = addfile("/work/home/L.xhm/LBM_standard/Fr/land.nc","r")
   land   = fi_lnd->land({50:75}, {150:200})

   Data_l = Data
   Data_o = Data
   do i = 0,dimsizes(Data(:,0,0,0,0))-1
   do j = 0,dimsizes(Data(0,:,0,0,0))-1
   do k = 0,dimsizes(Data(0,0,:,0,0))-1
      Data_l(i,j,k,:,:) = where(land.gt..4, Data(i,j,k,:,:),  Data@_FillValue)  ; land
      Data_o(i,j,k,:,:) = where(land.gt..4, Data@_FillValue,  Data(i,j,k,:,:))  ; ocean
   end do
   end do
   end do

   Data_a  = AreaMean(Data) 

   dim    := dimsizes(Data_a)
   data_rt = new((/3,dim(0),dim(1),dim(2)/), typeof(Data_a))
   data_rt(0,:,:,:)  = Data_a 
   data_rt(1,:,:,:)  = AreaMean(Data_l)
   data_rt(2,:,:,:)  = AreaMean(Data_o) 

   data_rt!3   = "lev" 
   data_rt&lev = Data&lev

   return(data_rt) 
end


undef("SigTest")
function SigTest(ctl, def, var, yyB, yyE)
local var_ctl,var_def,dim,dim_name
begin
   var_ctl = ReadData(ctl, var, yyB, yyE)
   var_def = ReadData(def, var, yyB, yyE)

   var_dif     = dim_avg_n(var_def,1) - dim_avg_n(var_ctl,1)
   var_ctl_avg = dim_avg_n(var_ctl,1)

   var_prob = ttest(dim_avg_n(var_def,1), dim_variance_n(var_def,1), dimsizes(var_def(0,:,0,0)),\
                    dim_avg_n(var_ctl,1), dim_variance_n(var_ctl,1), dimsizes(var_ctl(0,:,0,0)), True, False)
   var_prob@_FillValue = default_fillvalue("float")

   var_prob1  = var_prob
   var_prob2  = var_prob
   var_prob1  = where(var_prob.lt.0.05, var_dif, var_prob@_FillValue)
   var_prob2  = where(var_prob.lt.0.1,  var_dif, var_prob@_FillValue)

   rt_data = (/var_ctl_avg, var_dif, var_prob1, var_prob2/)
   copy_VarCoords(var_ctl(0,0,:,:), rt_data(0,0,:,:))
   return(rt_data)
end


undef("PLOT")
procedure PLOT(ctl, def, yyB, yyE, name)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = name ;"Pic/"+name

   ff_name = name+".nc"
   if(.not.fileexists(ff_name))then
      heat  = SigTest(ctl, def, "T", yyB, yyE)
      vapor = SigTest(ctl, def, "Q", yyB, yyE)

      ff_op = addfile(ff_name, "c")
      ff_op->heat  = heat
      ff_op->vapor = vapor
   else
      ff_in = addfile(ff_name, "r")
      heat  = ff_in->heat
      vapor = ff_in->vapor
   end if


; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

;  mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = "DJFM"

   lev = heat&lev
 lev_cc = ispan(1000,100,100) 
 lev_ID = lev_cc
 do i = 0,dimsizes(lev_cc)-1
    lev_ID(i) = ind(lev_cc(i) .eq. lev)
 end do

  res  = PLOT_xy_res(.4, .6, "T (K)",  "Height (hPa)", 3., "black", 0)
  res@gsnDraw              = True
  res@tmXTOn               = False
  res@tmYROn               = False
  res@tmYLMode             = "Explicit"
  res@tmYLValues           = log(lev(lev_ID))
  res@tmYLLabels           = lev(lev_ID)

  res@trYReverse           = True
  res@trYMinF              = log(100.)
  res@trYMaxF              = log(1000.)  
  res@trXMinF              = -.3
  res@trXMaxF              = .3 

  res@tiMainPosition       = "Left"
  res@tiXAxisFontColor     = "red"
  res@tmXBMajorLineColor   = "red"
  res@tmXBMinorLineColor   = "red"
  res@tmXBLabelFontColor   = "red"
  res@xyLineColors         = "red"
  res@xyMarkerColor        = "red"                    ; Marker color

  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Top"               ; Change location of
  res@pmLegendParallelPosF   = .768                  ; move units right
  res@pmLegendOrthogonalPosF = -0.26                ; move units down
  res@pmLegendWidthF         = 0.1                ; Change width and
  res@pmLegendHeightF        = 0.12                 ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .012       

  res@xyLineColors           = "red" 
  res@xyExplicitLegendLabels = "   T"

  res2 = res
  res2@tmXTOn              = True
  res2@tmXUseBottom        = False
  res2@tiXAxisString       = "Q (g/kg)"
  res2@tiXAxisFontColor    = "blue"
  res2@tmXTMajorLineColor  = "blue"
  res2@tmXTMinorLineColor  = "blue"
  res2@tmXTLabelFontColor  = "blue"
  res2@xyLineColors        = "blue"
  res2@xyMarkerColor       = "blue"                    ; Marker color

  res2@trXMinF             = -.01
  res2@trXMaxF             = .04 

  res2@pmLegendOrthogonalPosF = -0.5                ; move units down
  res2@xyLineColors           = "blue" 
  res2@xyExplicitLegendLabels = "   Q"

  res@tiMainString     = "(e) Vertical T & Q" 
 plot2 = gsn_csm_x2y(wks, heat(1,0,0,:), vapor(1,0,0,:)*1000., log(lev), res, res2)

  res@pmLegendDisplayMode  = "NoCreate"
  res@xyMarkLineModes      =  "Markers"
  res@xyMarkers            =  14                      ; choose type of marker
  res@xyMarkerSizeF        = 0.01

  res2@pmLegendDisplayMode = res@pmLegendDisplayMode
  res2@xyMarkLineModes     = res@xyMarkLineModes
  res2@xyMarkers           = res@xyMarkers
  res2@xyMarkerSizeF       = res@xyMarkerSizeF
 plot2 = gsn_csm_x2y(wks, heat(2,0,0,:), vapor(2,0,0,:)*1000., log(lev), res, res2)

frame(wks)


;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end

begin
; procedure PLOT(ctl, def, var, MinLat,MaxLat,MinLon,MaxLon, tiXString, yyB, yyE, name)

   ctl = "B2000_CAM5_CTR_"
   yyB = "0011"
   yyE = "0030"

   def = "B2000_CAM5_MC2BS_"

PLOT(ctl, def, yyB,yyE,  "Fig4e")

end
