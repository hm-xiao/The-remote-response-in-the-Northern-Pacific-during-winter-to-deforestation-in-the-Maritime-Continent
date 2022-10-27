load "Func/plot_map_res.ncl"
load "Func/above_terrain.ncl"

undef ("ReadData")
function ReadData(case, fi, Var, lev, yyB, yyE)
local Data,data,data_dt,data_avg,dim,dim_name
begin
   pwd = "/work/L.xhm/anal_ensemble/Comp/Data/"
   yyBI = toint(yyB)
   yyEI = toint(yyE)

   fis = systemfunc("ls "+pwd+case+"*/"+fi+"_0001_00?0.nc")
   fin = addfiles(fis, "r")
   ListSetType (fin, "join")

   if (Var.eq."precip")then
      PRECC = fin[:]->PRECC(:, yyBI*12-12:yyEI*12-1, : ,:)  
      PRECL = fin[:]->PRECL(:, yyBI*12-12:yyEI*12-1, : ,:)  
      Data  = (PRECC+PRECL)*30.*24.*3600.*1000.
      copy_VarCoords(PRECC, Data)
   else if( any(Var.eq.(/"U","V","VT","VQ","OMEGA","Z3","T","UVEL","VVEL","WVEL"/)) )then
      Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1, {lev}, : ,:)
   else
      Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1, : ,:)
   end if
   end if

   delete([/fis, fin/])
   if(any(Var.eq.(/"PSL","PS"/)))then
      data = Data/100.
      copy_VarCoords(Data,data)
      Data = data
      delete(data)
   end if

   if(Var.eq."Z3")then
      data = dim_rmvmean_n_Wrap(Data, 3)
      Data = data
      delete(data)
   end if

   dim  = dimsizes(Data)
   dim_name = getvardims(Data)

   Data@_FillValue     = default_fillvalue("float")
   Data@missing_value  = default_fillvalue("float")
   Data                = AboveTerrain(case, yyB, yyE, Data, lev)
;   Data(:,0,:,:) = default_fillvalue("float")
;   Data(:,dim(1)-2:,:,:) = default_fillvalue("float")

   data = reshape(Data, (/dim(0)*dim(1), dim(2), dim(3)/))
   data = runave_n(data, 4, 0, 0)

   res_data  = reshape(data(6-1:dim(0)*dim(1)-7-1,:,:), (/dim(0)*dim(1)/12-1, 12, dim(2),dim(3)/)) ;  #####
   res_data!2 = dim_name(2)
   res_data!3 = dim_name(3)
   res_data&$dim_name(2)$ = Data&$dim_name(2)$
   res_data&$dim_name(3)$ = Data&$dim_name(3)$
   return( res_data )
end

undef("SigTest")
function SigTest(ctl, def, fi, var, lev, yyB, yyE)
local var_ctl,var_def,dim,dim_name 
begin
   var_ctl = ReadData(ctl, fi, var, lev, yyB, yyE)
   var_def = ReadData(def, fi, var, lev, yyB, yyE)

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


undef("PLOT")
procedure PLOT(ctl, def,  ff1,var1, ff2,var2, ff3,var3, lev,  Levels,NN,RefMagnitudeF, LeftString, name, yyB, yyE)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = "Pic/"+name

   data1 = SigTest(ctl, def, ff1, var1, lev, yyB, yyE)
   data2 = SigTest(ctl, def, ff2, var2, lev, yyB, yyE)
   data3 = SigTest(ctl, def, ff3, var3, lev, yyB, yyE)

plot_data1 = data1[1]
plot_data2 = data2[1]
plot_data3 = data3[1]

lat = plot_data1&lat
lon = plot_data1&lon

; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  ;mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = (/"JJA", "SON", "DJF", "MAM"/)
  plot  = new(dimsizes(mm), graphic)
  plotT = new(dimsizes(plot_data1), graphic)
  plotT2= new(dimsizes(plot_data1), graphic)
  dum1  = new(dimsizes(mm), graphic)

  res = PLOT_map_res(Levels, NN, -10., 75., 80., 240.) 

  resline = PLOT_boxline_res(3., "purple")  ;LineThicknessF, LineColor

  resvc = PLOT_vc_res( 1.5, 0.018, 0.025, RefMagnitudeF, .08)  
  resvc@vcMonoFillArrowFillColor = False
  resvc@vcLevelSelectionMode     =  "ExplicitLevels"

  VcLevels  = ispan(1,100,1)/100.
  resvc@vcLevels                 = VcLevels
  LevelPlaette  = new(dimsizes(VcLevels), string)
  LevelPlaette  = "black"
  LevelPlaette(ind(VcLevels.lt..2)) = "white"
  resvc@vcLevelPalette           = LevelPlaette


  resdot  = PLOT_sig_res("grey66", .007)
  resdot2 = PLOT_sig_res("grey66", .003)

    do i = 0,0; dimsizes(mm)-1
       res@gsnLeftString   = LeftString
       ii = i+7
       plot(i)   = gsn_csm_contour_map_ce(wks, plot_data1(ii,:,:), res)
       dum1(i)   = gsn_add_polyline(wks, plot(i),(/90,90,150,150,90/),(/-10,10,10,-10,-10/), resline)
       aa        = gsn_csm_vector(  wks, plot_data2(ii,:,:),  plot_data3(ii,:,:), resvc)
       overlay(plot(i),  aa)
       delete([/aa/])

       do j = 0,dimsizes(data1[2](ii,:,0))-1
          do k = 0,dimsizes(data1[2](ii,0,:))-1
             if(.not.ismissing(data1[2](ii,j,k)))then
                plotT(i,j,k)  = gsn_add_polymarker(wks, plot(i), lon(k), lat(j), resdot)
             end if
             if(.not.ismissing(data1[3](ii,j,k)))then
                plotT2(i,j,k) = gsn_add_polymarker(wks, plot(i), lon(k), lat(j), resdot2)
             end if
          end do
       end do

    end do

  resP = True
;  resP@gsnPanelMainFontHeightF = 0.01
;  resP@gsnPanelLabelBar     = True
;  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller


;  gsn_panel(wks, plot, (/2,2/), resP)
  gsn_panel(wks, plot(0), (/1,1/), resP)

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

Levels_TS =  (/-3.6, -2.8, -2., -1.6, -1.2, -.8, -.4, -.2,\
                 .2,   .4,  .8,  1.2,  1.6,  2., 2.8, 3.6/)/4. 

Levels_PS =  (/-.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2,\
                .2,  .3,  .4,  .5,  .6,  .7,  .8,  .9/)
           
Levels_Z3 = (/-5., -4., -3., -2.5, -2., -1.5, -1., -.5,\
               .5,  1., 1.5,   2., 2.5,   3.,  4.,  5./) *2.

Levels_W  = (/-1.2,- .8, -.6,-.5,-.4,-.3, -.2,-.1,\
                .1,  .2,  .3, .4, .5, .6, .8, 1.2/) /100.

Levels_lh  = (/-1.6, -1.2, -1., -.8, -.6, -.4, -.2, -.1,\
                 .1,   .2,  .4,  .6,  .8,  1.,  1.2,  1.6/)
          
   do i = 0,1
      folder = where(i.eq.0, "BS", "C4")
      def    = where(i.eq.0, "B2000_CAM5_MC2BS_", "B2000_CAM5_MC2C4_")

PLOT(ctl, def,   "Z3","Z3", "U","U", "V","V", "200", Levels_Z3,8, 1., "(b) 200hPa Geopotential Height (m), u-v (m/s)", folder+"/200z3_uv_map_cli_dif", yyB, yyE)
exit

PLOT(ctl, def, "cam","precip", "U","U", "V","V","1000",  Levels_lh*10.,8,.5, "Shading: precip(mm/mon)", folder+"/precip_uv_map_cli_dif", yyB, yyE)

exit
;PLOT(ctl, def,   "Z3","Z3", "U","U", "V","V", "200", Levels_Z3,8, 2., "Shading: Z3(m); vector:200hPa u&v", folder+"/200z3_uv_map_cli_dif", yyB, yyE)

;PLOT(ctl, def,  "OMEGA","OMEGA", "U","U", "V","V", "200",Levels_W,8, 2.,  "Shading: 200hPa OMEGA(Pa/s); vector:200hPa u&v", folder+"/200w_uv_map_cli_dif", yyB, yyE)
;PLOT(ctl, def,  "OMEGA","OMEGA", "U","U", "V","V", "500",Levels_W,8, 2.,  "Shading: 500hPa OMEGA(Pa/s); vector:500hPa u&v", folder+"/500w_uv_map_cli_dif", yyB, yyE)
;PLOT(ctl, def,  "OMEGA","OMEGA", "U","U", "V","V", "850",Levels_W,8, 2.,  "Shading: 850hPa OMEGA(Pa/s); vector:850hPa u&v", folder+"/850w_uv_map_cli_dif", yyB, yyE)

delete(def)
   end do

end
