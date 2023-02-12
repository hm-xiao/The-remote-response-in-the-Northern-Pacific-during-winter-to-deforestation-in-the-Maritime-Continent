load "plot_map_res.ncl"

undef("PLOT")
procedure PLOT(area, bs_file, LeftString, var1,var2, FigName)
begin
   name = area+"_t42l20_"+bs_file
   
ff_name = FigName+".nc"
if(.not.fileexists(ff_name))then

    fin = "/work/home/L.xhm/LBM_standard/Out/"+area+"/"+name
    ff  = addfile(fin+".nc", "r")
     p  = rm_single_dims(ff->p)
     t  = rm_single_dims(ff->t(:,1,:,:))

   ff_op = addfile(ff_name, "c")
   ff_op->p=p
   ff_op->t=t
else
   ff_in = addfile(ff_name, "r")
   p = ff_in->p
   t = ff_in->t
end if


;=============================
; PLOT
;=============================
   figtype = "oldps"; "eps"
   figname = FigName 
    wks = gsn_open_wks(figtype,figname)


Levels_TS = (/-.5, -.4, -.3, -.25, -.2, -.15,  -.1, -.05,\
              .05,  .1, .15,   .2, .25,   .3,   .4, .5/)
Levels_slp  = (/-2.,-1.6, -1.2, -1., -.8, -.6, -.4, -.2, -.1,\
                .1,   .2,   .4,  .6,  .8,  1., 1.2, 1.6,  2./)

   res   = PLOT_map_res(Levels_TS, 8,  -10., 75., 80., 240. ) 
   resct = PLOT_contour_res(Levels_slp,  "black")

   res@gsnLeftString        = "(b) LBM 950hPa Air Temperature (K) & Surface Pressure (hPa), Avg: Day 30-59"
    plot = gsn_csm_contour_map_ce(wks, dim_avg_n_Wrap(t(20:49,:,:),0), res)
      aa = gsn_csm_contour(wks, dim_avg_n_Wrap(p(20:49,:,:),0), resct)
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
fr = (/"MC"/)
bs = (/"cli"/)

PLOT(fr,  bs,  "", "p","t","Fig4b")

end
     
