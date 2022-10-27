load "/work/L.xhm/anal_ensemble/Func/plot_map_res.ncl"

undef("PLOT")
procedure PLOT(area, bs_file, LeftString, var1,var2, FigName)
begin
   name = area+bs_file+"_t42l20_60dd"
   system("mkdir -p Pic/"+name) ;+"/png/")
   
   ntim = 59
   nvar = 7
   nlev = 20
   nlat = 64
   nlon = 128

   fin = "/work/home/L.xhm/LBM_standard/Out/"+area+"/"+name

  UNDEF = -999.                          ; UNDEF
  data  = new ( (/ntim,nvar,nlev,nlat,nlon/), float, UNDEF)
   slp  = new ( (/ntim,          nlat,nlon/), float, UNDEF)
  setfileoption("bin","ReadByteOrder","BigEndian")                                         

;  rec = fbinnumrec(fin+".grd")  8319 = 59*20*7 + 59
  k = 0
  do nt = 0,ntim-1                          ; read each record: store in x
  do nv = 0,nvar-1
  do nl = 0,nlev-1
     data(nt,nv,nl,:,:) = fbinrecread(fin+".grd", k, (/nlat,nlon/), "float")
     k = k+1
  end do
  end do
      slp(nt,:,:)       = fbinrecread(fin+".grd", k, (/nlat,nlon/), "float")
     k = k+1
  end do

  ff  = addfile(fin+".nc", "r")
;  p   = slp
;  t   = data(:,5,0,:,:)
;  copy_VarCoords(rm_single_dims(ff->p), p)
;  copy_VarCoords(rm_single_dims(ff->p), t)
   p  = rm_single_dims(ff->p)
   t  = rm_single_dims(ff->t(:,1,:,:))

;=============================
; POT
;=============================
   figtype = "oldps"; "eps"
   figname = "Pic/"+name+"/"+FigName
    wks = gsn_open_wks(figtype,figname)

Levels_TS = (/-.5, -.4, -.3, -.25, -.2, -.15,  -.1, -.05,\
              .05,  .1, .15,   .2, .25,   .3,   .4, .5/)
Levels_slp  = (/-2.,-1.6, -1.2, -1., -.8, -.6, -.4, -.2, -.1,\
                .1,   .2,   .4,  .6,  .8,  1., 1.2, 1.6,  2./)


  res = PLOT_map_res(Levels_TS, 8,  -10., 75., 80., 240.)
;  res@lbLabelBarOn       = False            ; turn off individual cb's


   resct = PLOT_contour_res(Levels_slp,  "black")

;;;;  Average PLOT------------------------------
    res@gsnLeftString        = "(b) LBM 950hPa Air Temperature (K) & Surface Pressure (hPa), Avg: Day 30-59"
       plot = gsn_csm_contour_map_ce(wks, dim_avg_n_Wrap(t(29:58,:,:),0), res)
       aa = gsn_csm_contour(wks, dim_avg_n_Wrap(p(29:58,:,:),0), resct)
       overlay(plot, aa)

  resP = True
  gsn_panel(wks, plot, (/1,1/), resP)


;========================Transform Image File========================
 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if
end

begin
;procedure PLOT(area, bs_file, LeftString, var1,var2, FigName)
event = (/"cli"/) ;, "nina", "nino"/)
mm = (/"Dec"/) ;(/"Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb","Mar", "Apr", "May"/)

do i = 0,dimsizes(event)-1
;do j = 0,dimsizes(area)-1
do k = 0,dimsizes(mm)-1
   PLOT("MC_CESM_lnd_DJFM",   "", "", "p","t","slp_temp")
end do
;end do
end do


end
     
