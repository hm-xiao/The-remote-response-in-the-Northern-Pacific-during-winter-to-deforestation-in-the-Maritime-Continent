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

   Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1,:,:); {50:75}, {140:190})

   if(any(Var.eq.(/"FLNS","LHFLX","SHFLX"/)))then
      data = -1.*Data
      copy_VarCoords(Data,data)
      Data = data
      delete(data)
   end if

   if(any(Var.eq.(/"CLDTOT","CLDLOW","CLDMED","CLDHGH"/)))then
      data = Data*10.^2
      copy_VarCoords(Data,data)
      Data = data
      delete(data)
   end if

   
   comp_data = Data ;(:,:,{-20:}, {60:260})

   dim  = dimsizes(comp_data)
   dim_name = getvardims(comp_data)

   comp_data@_FillValue     = default_fillvalue("float")
   comp_data@missing_value  = default_fillvalue("float")

   data = reshape(comp_data, (/dim(0)*dim(1), dim(2), dim(3)/))
   data = runave_n(data, 4, 0, 0)

   res_data  = reshape(data(6-1:dim(0)*dim(1)-7-1,:,:), (/dim(0)*dim(1)/12-1, 12, dim(2),dim(3)/)) ;  #####
   res_data!2 = dim_name(2)
   res_data!3 = dim_name(3)
   res_data&$dim_name(2)$ = comp_data&$dim_name(2)$
   res_data&$dim_name(3)$ = comp_data&$dim_name(3)$
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
procedure PLOT(ctl, def,  ff1,var1,lev,  Levels,NN, LeftString, name, yyB, yyE)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = name

   ff_name = name+".nc"
   if(.not.fileexists(ff_name))then
      data1 = SigTest(ctl, def, ff1, var1, lev, yyB, yyE)

      ff_op = addfile(ff_name, "c")
      ff_op->data1 = data1
   else
      ff_in = addfile(ff_name, "r")
      data1 = ff_in->data1
   end if

plot_data1 = data1(1,:,:,:)

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

  if(var1.eq."TS")then
     res = PLOT_map_res(Levels, NN,  -10., 10., 90., 150. )
  else
     res = PLOT_map_res(Levels, NN,   50., 75., 150., 200. )
     res@lbOrientation = "vertical"
  end if

  resline = PLOT_boxline_res(1., "purple")

  resdot  = PLOT_sig_res("grey66", .007)
  resdot2 = PLOT_sig_res("grey66", .003)

    do i = 0,dimsizes(mm)-1
       res@gsnLeftString   = LeftString
       ii = i+7
       plot(i)   = gsn_csm_contour_map_ce(wks, plot_data1(ii,:,:), res)
       dum1(i)  = gsn_add_polyline(wks, plot(i),(/90,90,150,150,90/),(/-10,10,10,-10,-10/), resline)

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

Levels_lh  = (/-2.4, -1.8, -1.2,  -.9, -.6, -.4,\
                 .4,   .6,   .9,  1.2, 1.8, 2.4/)
Levels_cloud = (/  -.1, -.08, -.06, -.04, -.02, -.01,\
                   .01,  .02,  .04,  .06,  .08,   .1/) *10.
Levels_TS2 = (/-5., -4., -3., -2.5, -2., -1.5, -1., -.5,\
                .5,  1., 1.5,   2., 2.5,   3.,  4., 5./)

      def    = "B2000_CAM5_MC2BS_"


PLOT(ctl, def, "cam","FLDS",  "sfc", Levels_lh,   6, " (c) Downwelling Long Wave Radiaton (W/m~S~2~N~)", "Fig4c", yyB, yyE)
PLOT(ctl, def, "cam","CLDHGH","sfc", Levels_cloud,6, " (d) Vertically Integrated High Cloud (%)",        "Fig4d", yyB, yyE)
PLOT(ctl, def, "cam","TS","    sfc",  Levels_TS2, 8, " (a) Skin temperature (K)",                        "FigS1a", yyB,yyE)

PLOT(ctl, def, "cam","FSDS",  "sfc", Levels_lh,   6, " (a) Downwelling Short Wave Radiaton (W/m~S~2~N~)","FigS5a", yyB, yyE)
PLOT(ctl, def, "cam","FSNS",  "sfc", Levels_lh,   6, " (b) Net Short Wave Radiaton (W/m~S~2~N~)",        "FigS5b", yyB, yyE)

end
