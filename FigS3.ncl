load "../Func/plot_map_res.ncl"
load "../Func/above_terrain.ncl"

undef ("ReadData")
function ReadData(case, fi, Var, lev, yyB, yyE)
local Data,data,data_dt,data_avg,dim,dim_name
begin
   pwd = "/work/L.xhm/anal_mean/Comp/Data/"
   yyBI = toint(yyB)
   yyEI = toint(yyE)

   fis = systemfunc("ls "+pwd+case+"*/"+fi+"_0001_00?0.nc")
   fin = addfiles(fis, "r")
   ListSetType (fin, "join")

   Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1, {lev}, : ,:)
;   Data = AboveTerrain(case, yyB, yyE, Data, lev)

   dim  = dimsizes(Data)
   dim_name = getvardims(Data)

   Data@_FillValue     = default_fillvalue("float")
   Data@missing_value  = default_fillvalue("float")
;   Data(:,0,:,:) = default_fillvalue("float")
;   Data(:,dim(1)-2:,:,:) = default_fillvalue("float")

   data = reshape(Data, (/dim(0)*dim(1), dim(2), dim(3)/))
   data = runave_n(data, 4, 0, 0)

   res_data  = reshape(data(6-1:dim(0)*dim(1)-7-1,:,:), (/dim(0)*dim(1)/12-1, 12, dim(2),dim(3)/)) ;  #####
   res_data := res_data(:,7:7,:,:)
   res_data!2 = dim_name(2)
   res_data!3 = dim_name(3)
   res_data&$dim_name(2)$ = Data&$dim_name(2)$
   res_data&$dim_name(3)$ = Data&$dim_name(3)$
   return( res_data )
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
      
   dim  = dimsizes(var_ctl)
   dim_name = getvardims(var_ctl)

   var_dif!1 = dim_name(2)
   var_dif!2 = dim_name(3)
   var_dif&$dim_name(2)$ = var_ctl&$dim_name(2)$
   var_dif&$dim_name(3)$ = var_ctl&$dim_name(3)$
   copy_VarCoords(var_dif,var_prob1)
   copy_VarCoords(var_dif,var_prob2)
   copy_VarCoords(var_dif,var_ctl_avg)
   return([/var_ctl_avg, var_dif, var_prob1, var_prob2/])
end
   
undef("CompData")
function CompData(case, lev, yyB, yyE)
begin
   uu  = ReadData(case, "U", "U", lev,  yyB, yyE)
   vv  = ReadData(case, "V", "V", lev,  yyB, yyE)
   sfvr = uv2sfvpG(uu, vv)
   sf = sfvr(0,:,:,:,:)*10^-6
   vp = sfvr(1,:,:,:,:)*10^-5
   copy_VarCoords(uu, sf)
   copy_VarCoords(uu, vp)

   div = uv2dvG_Wrap(uu, vv)
   ud  = new ( dimsizes(uu), typeof(uu), uu@_FillValue )
   vd  = new ( dimsizes(vv), typeof(vv), vv@_FillValue )

  dv2uvg(div,ud,vd) ; div  ==> divergent  wind components
  copy_VarCoords(uu, ud)
  copy_VarCoords(uu, vd)

  ud@long_name  = "Zonal Divergent Wind"
  ud@units      = uu@units
  vd@long_name  = "Meridional Divergent Wind"
  vd@units      = vv@units
   return([/uu, vv, ud, vd, sf, vp/])
end

undef("PLOT")
procedure PLOT(ctl, def, name, yyB, yyE, lev, RefMagnitudeF, Level, LeftString)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = "Pic/"+name

  data_ctl = CompData(ctl, lev, yyB, yyE)
  data_def = CompData(def, lev, yyB, yyE)

plot_u  = SigTest(data_ctl[2], data_def[2])
plot_v  = SigTest(data_ctl[3], data_def[3])
plot_vp = SigTest(dim_rmvmean_n_Wrap(data_ctl[5],3), dim_rmvmean_n_Wrap(data_def[5],3))
;plot_vp = SigTest(data_ctl[5], data_def[5])

lat = plot_vp[1]&lat
lon = plot_vp[1]&lon


; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  ;mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = (/"JJA", "SON", "DJF", "MAM"/)
  plot   = new((/dimsizes(mm)/), graphic)
  plotT1 = new(dimsizes(plot_u[2]), graphic)
  plotT2 = new(dimsizes(plot_u[2]), graphic)
  dum1  = new(dimsizes(mm), graphic)

  res = PLOT_map_res(Level, 9, -10, 75, 80, 240)

  resline = PLOT_boxline_res(3., "purple")  ;LineThicknessF, LineColor

    res2 = res
    delete([/res2@cnFillPalette/])
    res2@cnFillPalette       = ReadColorMap2("NCV_jet", 8)

  resvc = PLOT_vc_res( 1.6, 0.016, 0.025, RefMagnitudeF, .08)  
  resvc@vcMonoFillArrowFillColor = False
  resvc@vcLevelSelectionMode     =  "ExplicitLevels"

  VcLevels  = ispan(1,100,1)/100.
  resvc@vcLevels                 = VcLevels
  LevelPlaette  = new(dimsizes(VcLevels), string)
  LevelPlaette  = "black"
  LevelPlaette(ind(VcLevels.lt.RefMagnitudeF/5.)) = "white"
  resvc@vcLevelPalette           = LevelPlaette

  resdot  = PLOT_sig_res("grey66", .007)
  resdot2 = PLOT_sig_res("grey66", .003)

  resP = True
;  resP@gsnPanelLabelBar     = True
;  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller

    do i = 0,0 ;dimsizes(mm)-1
;       ii = i*3+1
       ii = i
       resvc@vcRefMagnitudeF = RefMagnitudeF

       res@gsnLeftString     = LeftString
       plot(i)  = gsn_csm_contour_map_ce(wks, plot_vp[1](ii,:,:), res)
       dum1(i)   = gsn_add_polyline(wks, plot(i),(/90,90,150,150,90/),(/-10,10,10,-10,-10/), resline)
       aa        = gsn_csm_vector(  wks, plot_u[1](ii,:,:),  plot_v[1](ii,:,:), resvc)
       overlay(plot(i),  aa)
       delete([/aa/])

       do j = 0,dimsizes(lat)-1
          do k = 0,dimsizes(lon)-1
             if(.not.ismissing(plot_vp[2](ii,j,k)))then
                plotT1(i,j,k) = gsn_add_polymarker(wks, plot(i), lon(k), lat(j), resdot)
             end if
             if(.not.ismissing(plot_vp[3](ii,j,k)))then
                plotT2(i,j,k) = gsn_add_polymarker(wks, plot(i), lon(k), lat(j), resdot2)
             end if
          end do
       end do

    end do

    gsn_panel(wks, (/plot/),  (/1,1/),  resP)

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

   do i = 0,0
      folder = where(i.eq.0, "BS", "C4")
      def    = where(i.eq.0, "B2000_CAM5_MC2BS_", "B2000_CAM5_MC2C4_")
      PLOT(ctl, def, folder+"/850vp_uv_map_cli_dif", yyB, yyE, 850, .25,\
(/-2.5, -2., -1.5, -1.2, -1., -.8, -.6, -.5, .5, .6, .8, 1., 1.2, 1.5, 2., 2.5/),"(a) 850hPa Velocity Potential (m~S~2~N~/s), u-v (m/s)")

      PLOT(ctl, def, folder+"/200vp_uv_map_cli_dif", yyB, yyE, 200, .5,\
(/-4.5, -4., -3.5, -3., -2.5, -2., -1.5, -1., 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5/),"(b) 200hPa Velocity Potential (m~S~2~N~/s), u-v (m/s)")

   end do
end
