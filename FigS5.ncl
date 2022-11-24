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

   if(Var.eq."Rnet")then
      FLNS = fin[:]->FLNS(:, yyBI*12-12:yyEI*12-1,:,:)
      Data = -FLNS + fin[:]->FSNS(:, yyBI*12-12:yyEI*12-1,:,:)
      copy_VarCoords(FLNS,Data)
   else if(Var.eq."FLUS")then
      FLNS = fin[:]->FLNS(:, yyBI*12-12:yyEI*12-1,:,:)
      Data =  -( FLNS + fin[:]->FLDS(:, yyBI*12-12:yyEI*12-1,:,:) )
      copy_VarCoords(FLNS,Data)
   else if(any(Var.eq.(/"T","VT","VQ"/)))then
      Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1, {lev}, :,:)

   else if(any(Var.eq.(/"precip","snow","snow_ratio"/)))then
      PRECC = fin[:]->PRECC(:, yyBI*12-12:yyEI*12-1, : ,:)
      PRECL = fin[:]->PRECL(:, yyBI*12-12:yyEI*12-1, : ,:)
      PRECT = PRECC+PRECL
      Data  = PRECT*30.*24.*3600.*1000.
      if(any(Var.eq.(/"snow","snow_ratio"/)))then
         PRECSC = fin[:]->PRECSC(:, yyBI*12-12:yyEI*12-1, : ,:)
         PRECSL = fin[:]->PRECSL(:, yyBI*12-12:yyEI*12-1, : ,:)
         PRECST = PRECSC+PRECSL
         Data   = PRECST*30.*24.*3600.*1000.
         if(Var.eq."snow_ratio")
            PRECT = where(PRECT.eq.0., default_fillvalue("float"), PRECT)
            Data  = PRECST/PRECT
         end if
      end if
      copy_VarCoords(PRECC, Data)

   else if(Var.eq."evapotr_flux")then   ; evapotr_mm: QSOIL+QVEGE+QVEGT
      FCEV = fin[:]->FCEV(:, yyBI*12-12:yyEI*12-1, : ,:)
      FCTR = fin[:]->FCTR(:, yyBI*12-12:yyEI*12-1, : ,:)
      FGEV = fin[:]->FGEV(:, yyBI*12-12:yyEI*12-1, : ,:)
      Data  = -(FCEV+FCTR+FGEV)    ; Down: Positive
      copy_VarCoords(FCEV, Data)

   else
      Data = fin[:]->$Var$(:, yyBI*12-12:yyEI*12-1,:,:); {50:75}, {140:190})
   end if
   end if
   end if
   end if
   end if

   if(Var.eq."TMQ")then
      data = Data/1000.*1000.    ; kg/m^2 -> mm
      copy_VarCoords(Data,data)
      Data = data
      delete(data)
   end if

   if(Var.eq."QREFHT" .or. Var.eq."VQ")then
      data = Data*10.^3
      copy_VarCoords(Data,data)
      Data = data
      delete(data)
   end if

   if(any(Var.eq.(/"FLNS","LHFLX","SHFLX"/)))then
      data = -Data
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
   comp_data = AboveTerrain(case, yyB, yyE, comp_data, lev)

   dim  = dimsizes(comp_data)
   dim_name = getvardims(comp_data)

   comp_data@_FillValue     = default_fillvalue("float")
   comp_data@missing_value  = default_fillvalue("float")
;   comp_data(:,0,:,:) = default_fillvalue("float")
;   comp_data(:,dim(1)-2:,:,:) = default_fillvalue("float")

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
procedure PLOT(ctl, def,  ff1,var1,lev,  Levels,NN, LeftString, name, yyB, yyE)
local data1,data2,data3
begin
   figtype = "oldps"
   figname = "Pic/"+name

   data1 = SigTest(ctl, def, ff1, var1, lev, yyB, yyE)

plot_data1 = data1[1]
;plot_data1 = where(data1[3].lt.0.1,                     data1[1], data1[1]@_FillValue)

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

  ;res = PLOT_map_res(Levels, NN,  -10., 75., 80., 240. )
  res = PLOT_map_res(Levels, NN,   50., 75., 150., 200. )
  res@lbLabelBarOn         = False
;  res@lbOrientation = "vertical"
;  res@lbLabelFontHeightF   = 0.016               ; make labels smaller

  resline = PLOT_boxline_res(1., "purple")

  resdot  = PLOT_sig_res("grey66", .007)
  resdot2 = PLOT_sig_res("grey66", .003)

    do i = 0,0 ;dimsizes(mm)-1
;       res@tiMainString     = "mon:"+mm(i)
       res@gsnLeftString   = LeftString
;       ii = i*3+1
       ii = i+7
       plot(i)   = gsn_csm_contour_map_ce(wks, plot_data1(ii,:,:), res)
       dum1(i)  = gsn_add_polyline(wks, plot(i),(/90,90,150,150,90/),(/-10,10,10,-10,-10/), resline)

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
Levels_TS = (/-.5, -.4, -.35, -.3, -.25, -.2, -.15, -.1,\
               .1, .15,   .2, .25,   .3, .35,   .4, .5/)



   do i = 0,1
      folder = where(i.eq.0, "BS", "C4")
      def    = where(i.eq.0, "B2000_CAM5_MC2BS_", "B2000_CAM5_MC2C4_")

PLOT(ctl, def, "cam","FSDS","sfc",  Levels_lh,6, "(a) Downwelling Short Wave Radiaton (W/m~S~2~N~)",  folder+"/FSDS_map_cli_dif", yyB, yyE)
exit
PLOT(ctl, def, "cam","FSNS","sfc",  Levels_lh,6, "(b) Net Short Wave Radiaton (W/m~S~2~N~)",  folder+"/FSNS_map_cli_dif", yyB, yyE)
exit
PLOT(ctl, def, "cam","FLDS","sfc",  Levels_lh,6, "(c) Downwelling Long Wave Radiaton (W/m~S~2~N~)",  folder+"/FLDS_map_cli_dif", yyB, yyE)
PLOT(ctl, def, "cam","CLDHGH","sfc",  Levels_cloud,6, "(d) Vertically Integrated High Cloud (%)",  folder+"/CLDHGH_map_cli_dif", yyB, yyE)
exit
;PLOT(ctl, def, "cam","CLDLOW","sfc",  Levels_cloud,8, "(d) Vertically Integrated Low Cloud (%)",  folder+"/CLDLOW_map_cli_dif", yyB, yyE)
PLOT(ctl, def, "cam","CLDTOT","sfc",  Levels_cloud,8, "(d) Vertically Integrated Total Cloud (%)",  folder+"/CLDTOT_map_cli_dif", yyB, yyE)
PLOT(ctl, def, "cam","CLDMED","sfc",  Levels_cloud,8, "(d) Vertically Integrated Mid-level Cloud (%)",  folder+"/CLDMED_map_cli_dif", yyB, yyE)
exit

exit
   end do

end
