load "plot_map_res.ncl"
load "compute.ncl"

undef ("ReadData")
function ReadData(case, yyB, yyE)
local Data,data,data_dt,data_avg,dim,dim_name
begin
   pwd = "/work/L.xhm/anal_mean/Comp/Data/"
   yyBI = toint(yyB)
   yyEI = toint(yyE)

   fis = systemfunc("ls "+pwd+case+"*/OMEGA_0001_00?0.nc")
   fin = addfiles(fis, "r")
   ListSetType (fin, "join")

   Data = fin[:]->OMEGA(:, yyBI*12-12:yyEI*12-1, :, {-10:10}, {90:150})
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
   land   = fi_lnd->land({-10:10} ,{90:150})

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
function SigTest(ctl, def, yyB, yyE)
local var_ctl,var_def,dim,dim_name
begin
   var_ctl = ReadData(ctl, yyB, yyE)
   var_def = ReadData(def, yyB, yyE)
   var_prob = ttest(dim_avg_n(var_def,1), dim_variance_n(var_def,1), dimsizes(var_def(0,:,0,0)),\
                    dim_avg_n(var_ctl,1), dim_variance_n(var_ctl,1), dimsizes(var_ctl(0,:,0,0)), True, False)
   var_prob@_FillValue = default_fillvalue("float")

   var_prob1  = var_prob
   var_prob2  = var_prob
   var_prob1  = where(var_prob.lt.0.05, var_prob, var_prob@_FillValue)
   var_prob2  = where(var_prob.lt.0.1,  var_prob, var_prob@_FillValue)

   var_dif     = dim_avg_n(var_def,1) - dim_avg_n(var_ctl,1)
   var_ctl_avg = dim_avg_n(var_ctl,1)
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
      heat = SigTest(ctl, def, yyB, yyE)

      ff_op = addfile(ff_name, "c")
      ff_op->heat = heat
   else
      ff_in = addfile(ff_name, "r")
      heat  = ff_in->heat
   end if


; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

;  mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = "DJFM"
  plot2 = new((/  dimsizes(mm)/), graphic)

   lev = heat&lev
 lev_cc = ispan(1000,100,100) 
 lev_ID = lev_cc
 do i = 0,dimsizes(lev_cc)-1
    lev_ID(i) = ind(lev_cc(i) .eq. lev)
 end do

  RES  = PLOT_xy_res(.4, .6, "OMEGA (Pa/s)",  "Height (hPa)", 3., "black", 0)
  res  = RES
  res@tmXTOn               = False
  res@tmYROn               = False
  res@tmYLMode             = "Explicit"
  res@tmYLValues           = log(lev(lev_ID))
  res@tmYLLabels           = lev(lev_ID)

  res@trYReverse           = True
  res@trYMinF              = log(100.)
  res@trYMaxF              = log(1000.)  

  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .4                  ; move units right
  res@pmLegendOrthogonalPosF = -0.26                ; move units down
  res@pmLegendWidthF         = 0.072                ; Change width and
  res@pmLegendHeightF        = 0.12                 ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .012                 ; label font height

  delete(res@xyLineColors)
  res@xyLineColors           = (/"black",  "brown",    "blue"/)
  res@xyExplicitLegendLabels = (/"  All", "  Land", "  Ocean"/)

  res2 = res

    do i = 0,0;dimsizes(mm)-1
       ii = i;+6
       res2@gsnLeftString     = "(b)" ;"mon:"+mm(i)
      plot2(i) = gsn_csm_xy(wks, (/heat(1,0,ii,:), heat(1,1,ii,:),heat(1,2,ii,:)/), log(lev), res2)
      add_marker(wks, plot2(i), "black",  ind(heat(2,0,ii,:).lt..05),   heat(1,0,ii,:), log(lev))
      add_marker(wks, plot2(i), "brown",  ind(heat(2,1,ii,:).lt..05),   heat(1,1,ii,:), log(lev))
      add_marker(wks, plot2(i), "blue",   ind(heat(2,2,ii,:).lt..05),   heat(1,2,ii,:), log(lev))
    end do

  resP = True
;  resP@gsnPanelLabelBar     = True
;  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller

  gsn_panel(wks, plot2(0),  (/1,1/), resP)

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

PLOT(ctl, def, yyB,yyE,  "FigS1b")

end
