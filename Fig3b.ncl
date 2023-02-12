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

   Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1, {lev}, : ,:)

   dim  = dimsizes(Data)
   dim_name = getvardims(Data)

   Data@_FillValue     = default_fillvalue("float")
   Data@missing_value  = default_fillvalue("float")
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
   
undef("CompData")
function CompData(ctl, def, lev, yyB, yyE)
begin
   u_c  = ReadData(ctl, "U", "U", lev,  yyB, yyE)
   v_c  = ReadData(ctl, "V", "V", lev,  yyB, yyE)
   sfvr_c = uv2sfvpG(u_c, v_c)
   sf_c = sfvr_c(0,:,:,:,:)*10^-6
   vp_c = sfvr_c(1,:,:,:,:)*10^-5
   copy_VarCoords(u_c,sf_c)
   copy_VarCoords(u_c,vp_c)

   u_d  = ReadData(def, "U", "U", lev,  yyB, yyE)
   v_d  = ReadData(def, "V", "V", lev,  yyB, yyE)
   sfvr_d = uv2sfvpG(u_d, v_d)
   sf_d = sfvr_d(0,:,:,:,:)*10^-6
   vp_d = sfvr_d(1,:,:,:,:)*10^-5
   copy_VarCoords(u_d,sf_d)
   copy_VarCoords(u_d,vp_d)
   return([/u_c, u_d, v_c, v_d, sf_c, sf_d, vp_c, vp_d/])
end

undef("PLOT")
procedure PLOT(ctl, def, name, yyB, yyE)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = name

   ff_name = name+".nc"
   if(.not.fileexists(ff_name))then

      u_c_o= ReadData(ctl, "U", "U", 200,  yyB, yyE)
      u_d  = ReadData(def, "U", "U", 200,  yyB, yyE)
      u_c  = dim_rmvmean_n_Wrap(u_c_o, 3)
      u_d  = dim_rmvmean_n_Wrap(u_d, 3)
      u_comp =  SigTest(u_c, u_d)
      u_c_o_avg = dim_avg_n_Wrap(u_c_o, 0)

      ff_op = addfile(ff_name, "c")
      ff_op->u_comp=u_comp
      ff_op->u_c_o_avg=u_c_o_avg
   else
      ff_in = addfile(ff_name, "r")
      u_comp    = ff_in->u_comp
      u_c_o_avg = ff_in->u_c_o_avg
   end if

   lat = u_comp&lat
   lon = u_comp&lon

; Plot—————————————
   wks = gsn_open_wks(figtype, figname)

  ;mm = (/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May"/)
  mm = "DJFM" 
  plot  = new(dimsizes(mm), graphic)
  plotT = new(dimsizes(u_comp(1,:,:,:)), graphic)
  plotT2= new(dimsizes(u_comp(1,:,:,:)), graphic)


 Levels  = (/-1.5, -1.2, -.8, -.6, -.5, -.4, -.3, \
               .3,   .4,  .5,  .6,  .8, 1.2, 1.5/)
            

  res   = PLOT_map_res(Levels, 7,  -10., 75., 80., 240.)
  resvc = PLOT_vc_res( 1.8, 0.012, 0.03, .1, .08)  

  resdot  = PLOT_sig_res("grey66", .007)
  resdot2 = PLOT_sig_res("grey66", .003)

  res2 = PLOT_contour_res(ispan(-10,60,10), "black")
  res2@gsnContourLineThicknessesScale = 1.36

  resP = True
;  resP@gsnPanelLabelBar     = True
;  resP@lbLabelFontHeightF   = 0.012               ; make labels smaller

    do i  = 0,dimsizes(mm)-1
       ii = i+7

      res@gsnLeftString     = "(b) 200hPa Zonal Wind (m/s)"
       plot(i)  = gsn_csm_contour_map_ce(wks, u_comp(1,ii,:,:), res)
       cli      = gsn_csm_contour( wks, u_c_o_avg(ii,:,:),  res2)

       overlay(plot(i),  cli)

       do j = 0,dimsizes(lat)-1
          do k = 0,dimsizes(lon)-1
             if(.not.ismissing(u_comp(2,ii,j,k)))then
                plotT(i,j,k)  = gsn_add_polymarker(wks, plot(i), lon(k), lat(j), resdot)
             end if
          end do
       end do

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
PLOT(ctl, def, "Fig3b", yyB, yyE)

end
