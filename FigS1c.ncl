;  https://renqlsysu.github.io/2019/02/01/apparent_heat_source/
;  https://doi.org/10.1175/2007MWR2145.1    ,Fig.6
load "plot_map_res.ncl"
load "compute.ncl"

undef ("ReadData")
function ReadData(case, var, yyB, yyE)
local Data,data,data_dt,data_avg,dim,dim_name
begin
   pwd = "/work/L.xhm/anal_mean/Comp/Data/"

   fin = addfile(pwd+case+"apparent_heat_source_DJFM.nc","r")
   Data = fin->$var$(:, :, :, {-10:10} ,{90:150})

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

   Data_a  = AreaMean(Data)   *24.*60.*60.   ; K/s -> C/day

   dim     = dimsizes(Data_a)
   data_rt = new((/3,dim(0),dim(1),dim(2)/), typeof(Data_a))
   data_rt(0,:,:,:)  = Data_a 
   data_rt(1,:,:,:)  = AreaMean(Data_l) *24.*60.*60.   ; K/s -> C/day
   data_rt(2,:,:,:)  = AreaMean(Data_o) *24.*60.*60.   ; K/s -> C/day

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
      heat = SigTest(ctl, def, "heat", yyB, yyE)

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

  RES  = PLOT_xy_res(.4, .6, "Diabatic heating (K/Day)", "Height (hPa)", 3., "black", 0)
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
  res@pmLegendParallelPosF   = .76                  ; move units right
  res@pmLegendOrthogonalPosF = -0.26                ; move units down
  res@pmLegendWidthF         = 0.072                ; Change width and
  res@pmLegendHeightF        = 0.12                 ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .012                 ; label font height

  delete(res@xyLineColors)
  res@xyLineColors           = (/"green", "blue", "red", "purple", "black"/)
  res@xyExplicitLegendLabels = (/"  Longwave", "  Shortwave", "  T (moist)", "  T (diffusion)", "  Total Heat"/)  

  res2 = res
  delete([/res2@xyLineColors, res2@xyExplicitLegendLabels/])
  res2@xyLineColors          = (/"black",  "brown",    "blue"/)
  res2@xyExplicitLegendLabels= (/"  All", "  Land", "  Ocean"/)

    do i = 0,0;dimsizes(mm)-1
       ii = i;+6
       res2@gsnLeftString     = "(c)" ;"mon:"+mm(i)
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
;procedure PLOT(ctl, def, yyB, yyE, area, name)

   ctl = "B2000_CAM5_CTR_"
   yyB = "0011"
   yyE = "0030"

   def = "B2000_CAM5_MC2BS_"

PLOT(ctl, def, yyB,yyE,  "FigS1c")


end
