load "plot_map_res.ncl"

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

   if( any(Var.eq.(/"U","V","VT","VQ","OMEGA","Z3","T","UVEL","VVEL","WVEL"/)) )then
      Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1, {lev}, : ,:)
   else
      Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1, : ,:)
   end if

   if(any(Var.eq.(/"PSL","PS"/)))then
      data = Data/100.
      copy_VarCoords(Data,data)
      Data = data
      delete(data)
   end if

   dim  = dimsizes(Data)
   dim_name = getvardims(Data)

   Data@_FillValue     = default_fillvalue("float")
   Data@missing_value  = default_fillvalue("float")

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
   var_ctl_avg = dim_avg_n(var_ctl, 0)

   rt_data = (/var_ctl_avg, var_dif, var_prob1, var_prob2/)
   copy_VarCoords(var_ctl(0,0,:,:), rt_data(0,0,:,:))
   return(rt_data)
end
   

undef("PLOT")
procedure PLOT(ctl, def,  ff1,var1,lev1, ff2,var2, ff3,var3, lev, Levels, NN, cnlev, RefMagnitudeF, LeftString, name, yyB, yyE)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = name

   ff_name = name+".nc"
   if(.not.fileexists(ff_name))then
      data1 = SigTest(ctl, def, ff1, var1, lev1, yyB, yyE)
      data2 = SigTest(ctl, def, ff2, var2, lev,  yyB, yyE)
      data3 = SigTest(ctl, def, ff3, var3, lev,  yyB, yyE)

      ff_op = addfile(ff_name, "c")
      ff_op->data1 = data1
      ff_op->data2 = data2
      ff_op->data3 = data3
   else
      ff_in = addfile(ff_name, "r")
      data1 = ff_in->data1
      data2 = ff_in->data2
      data3 = ff_in->data3
   end if

plot_data1 = data1(1,:,:,:)
plot_data2 = data2(1,:,:,:)
plot_data3 = data3(1,:,:,:)
data_c_avg = data1(0,:,:,:)

lat = plot_data1&lat
lon = plot_data1&lon

; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  ;mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = "DJFM"
  plot  = new(dimsizes(mm), graphic)
  plotT = new(dimsizes(plot_data1), graphic)
  plotT2= new(dimsizes(plot_data1), graphic)
  dum1  = new(dimsizes(mm), graphic)

  res  = PLOT_map_res(Levels, NN,  -10., 75., 80., 240. ) 

  RES_ct = PLOT_contour_res(cnlev, "grey66")  

  resline = PLOT_boxline_res(3., "purple")

  resvc = PLOT_vc_res( 1.8, 0.02, 0.03, RefMagnitudeF, .08)
  resvc@vcMonoFillArrowFillColor = False
  resvc@vcLevelSelectionMode     =  "ExplicitLevels"

  VcLevels  = ispan(5,50,1)/100.
  resvc@vcLevels                 = VcLevels 
  LevelPlaette  = new(dimsizes(VcLevels), string) 
  LevelPlaette  = "black"
  LevelPlaette(ind(VcLevels.lt..1)) = "white"
  resvc@vcLevelPalette           = LevelPlaette

  resdot  = PLOT_sig_res("grey66", .007)
  resdot2 = PLOT_sig_res("grey66", .003)

    do i = 0,dimsizes(mm)-1
       res@gsnLeftString     = LeftString 
       ii = i+7
       plot(i)   = gsn_csm_contour_map_ce(wks, plot_data1(ii,:,:), res)
       dum1(i)   = gsn_add_polyline(wks, plot(i),(/90,90,150,150,90/),(/-10,10,10,-10,-10/), resline)

       cli       = gsn_csm_contour( wks, data_c_avg(ii,:,:),  RES_ct)
       aa        = gsn_csm_vector(  wks, plot_data2(ii,:,:),  plot_data3(ii,:,:), resvc)
;       overlay(plot(i),  cli)
       overlay(plot(i),  aa)
       delete([/cli, aa/])

       do j = 0,dimsizes(data1(2,ii,:,0))-1
          do k = 0,dimsizes(data1(2,ii,0,:))-1
             if(.not.ismissing(data1(2,ii,j,k)))then
                plotT(i,j,k)  = gsn_add_polymarker(wks, plot(i), lon(k), lat(j), resdot)
             end if
          end do
       end do

    end do

  resP = True
;  resP@gsnPanelLabelBar     = True
;  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller


;  gsn_panel(wks, plot, (/3,4/), resP)
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

Levels_PS =  (/-1., -.8, -.7, -.6, -.5, -.4, -.3, -.2,\
                .2,  .3,  .4,  .5,  .6,  .7,  .8,  1./)

   def = "B2000_CAM5_MC2BS_"

PLOT(ctl, def,  "cam","PSL", "sfc", "U","U", "V","V", 1000, Levels_PS, 8, fspan(980,1030,11), .5, "(a) Sea Level Pressure (hPa), 1000hPa u-v (m/s)", "Fig1a", yyB, yyE)

end
