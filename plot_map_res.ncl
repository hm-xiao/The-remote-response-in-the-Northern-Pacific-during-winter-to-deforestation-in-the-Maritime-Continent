undef("ReadColorMap2")
function ReadColorMap2(cmap, NN)
begin
   nn = NN*2+1
   icol = span_color_rgba(cmap, nn)
   icol(NN,:) = (/1,1,1,1/)
   return(icol)
end


undef("PLOT_map_res")
function PLOT_map_res(Levels, NN, MinLatF, MaxLatF, MinLonF, MaxLonF)
begin
    res                      = True             ; plot mods desired
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnAddCyclic         = False
    res@gsnLeftString        = ""
    res@tiMainFontHeightF    = 0.02

    res@cnFillOn             = True               ; turn on color for contours
;    res@cnFillMode           = "CellFill"
    res@cnLinesOn            = False              ; turn off contour lines
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    res@cnInfoLabelOn        = False

    res@cnLevelSelectionMode = "ExplicitLevels" 
    res@cnLevels             = Levels
    res@cnFillPalette        = ReadColorMap2("NCV_jet", NN)

;    res@lbLabelBarOn         = False
    res@lbLabelFontHeightF   = 0.01               ; make labels smaller

    res@tmXTOn               = False
    res@tmYROn               = False
    res@tmBorderThicknessF   = .8
    res@tmXBMajorThicknessF  = 1.2
    res@tmYLMajorThicknessF  = 1.2
    res@tmXBLabelFontHeightF = 0.012
    res@tmYLLabelFontHeightF = 0.012

    res@mpMinLatF            = MinLatF 
    res@mpMaxLatF            = MaxLatF
    res@mpMinLonF            = MinLonF
    res@mpMaxLonF            = MaxLonF
    res@mpCenterLonF         = (res@mpMinLonF + res@mpMaxLonF)/2.
    res@mpGeophysicalLineThicknessF = .8
    res@mpLandFillColor      = "white"
return(res)
end


undef("PLOT_contour_res")
function PLOT_contour_res(cnlev, color)
begin
    res              = True                         ; plot mods desired
    res@gsnDraw      = False                        ; don't draw yet
    res@gsnFrame     = False                        ; don't advance frame yet
    res@gsnLeftString= ""
    res@gsnContourLineThicknessesScale = 1.18
    res@gsnContourZeroLineThicknessF = 1.     ; doubles thickness of zero contour
    res@gsnContourNegLineDashPattern = 11
    res@gsnContourPosLineDashPattern = 0
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels             = cnlev
    res@cnLineLabelBackgroundColor = -1
    res@cnLineColor          = color 
    res@cnInfoLabelOn        = False
    res@cnInfoLabelFontColor = color 
    res@cnLineLabelFontColor = color
return(res)
end


undef("PLOT_boxline_res")
function PLOT_boxline_res(LineThicknessF, LineColor)
begin
    resline                  = True
    resline@gsLineThicknessF = LineThicknessF
    resline@gsLineColor      = LineColor
return(resline)
end


undef("PLOT_vc_res")
function PLOT_vc_res(FillArrowEdgeThicknessF, MinDistanceF, RefLengthF, RefMagnitudeF, RefAnnoOrthogonalPosF)
begin
    resvc                    = True
    resvc@gsnDraw            = False
    resvc@gsnFrame           = False
    resvc@gsnAddCyclic       = False    ; data already has cyclic point
    resvc@gsnLeftString      = ""
    resvc@gsnRightString     = ""

    resvc@vcGlyphStyle       = "FillArrow" ;"LineArrow"
    resvc@vcFillArrowEdgeThicknessF = FillArrowEdgeThicknessF
    resvc@vcFillArrowHeadXF         = .45
    resvc@vcFillArrowHeadYF         = .18

;    resvc@vcLineArrowThicknessF = LineArrowThicknessF
    resvc@vcMapDirection     = False
    resvc@vcMinDistanceF     = MinDistanceF
    resvc@vcRefLengthF       = RefLengthF 
    resvc@vcRefMagnitudeF    = RefMagnitudeF               ; define vector ref mag
    resvc@vcRefAnnoOrthogonalPosF = RefAnnoOrthogonalPosF 
    resvc@vcRefAnnoFontHeightF    = .008
return(resvc)
end


undef("PLOT_sig_res")
function PLOT_sig_res(MarkerColor, MarkerSizeF)
begin
    resdot                  = True
    resdot@gsMarkerColor    = MarkerColor
    resdot@gsMarkerSizeF    = MarkerSizeF 
    resdot@gsMarkerIndex    = 1
return(resdot)
end


undef("PLOT_xy_res")
function PLOT_xy_res(WidthF, HeightF, XAxisString, YAxisString, LineThicknessF, LineColors, DashPatterns) 
begin
    res                      = True              ; plot mods desired
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@vpWidthF             = WidthF            ; Change the aspect ratio, but 
    res@vpHeightF            = HeightF           ; make plot as large as possible.
    res@tiXAxisString        = XAxisString
    res@tiYAxisString        = YAxisString 
    res@xyLineThicknessF     = LineThicknessF
    res@xyLineColors         = LineColors
    res@xyDashPatterns       = DashPatterns

    res@tiYAxisFontColor     = LineColors
    res@tmYLMajorLineColor   = LineColors
    res@tmYLMinorLineColor   = LineColors
    res@tmYLLabelFontColor   = LineColors
return(res)
end


undef("add_marker")
procedure add_marker(wks,plot,MarkerColor,prob_judge,xx,yy)
local dumstr
begin
  resD                  = True
  resD@gsMarkerIndex    = 14
  resD@gsMarkerSizeF    = 0.01
  resD@gsMarkerColor    = MarkerColor

  dumstr               = unique_string("marker")
  if(any(.not.ismissing(prob_judge)))then
     plot@$dumstr$ = gsn_add_polymarker(wks, plot, xx(prob_judge), yy(prob_judge), resD)
  end if
end


undef("add_gon")
procedure add_gon(wks,plot,cname,oval,xx,yy)
local gsres, dumstr
begin
  gsres                = True
  gsres@gsFillColor    = cname
  gsres@gsFillOpacityF = oval
  dumstr               = unique_string("gon")
  plot@$dumstr$ = gsn_add_polygon(wks,plot,xx,yy,gsres)
end
