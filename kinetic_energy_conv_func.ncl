; Perturbation kinetic energy conversion Reference
; https://link.springer.com/article/10.1007/s00382-022-06281-6

undef ("KPE")
function KPE(u_bar, v_bar, u_pri, v_pri)
begin
  uu_pri_avg = u_pri^2
  vv_pri_avg = v_pri^2
  uv_pri_avg = u_pri*v_pri

   lat = tofloat(u_bar&lat)
   lon = tofloat(u_bar&lon)

dU_bardx = u_bar
dV_bardx = u_bar
   dlon = (lon(2)-lon(1))*0.0174533
  do nl = 0,dimsizes(lat)-1                      ; loop over each latitude
     dx = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
     dU_bardx(nl,:)  = center_finite_diff_n(u_bar(nl:nl,:), dx, False,0,1) ;nl:nl is to preserve the 3D structure
     dV_bardx(nl,:)  = center_finite_diff_n(v_bar(nl:nl,:), dx, False,0,1) ;nl:nl is to preserve the 3D structure
  end do

     dU_bardy   = center_finite_diff_n(u_bar, lat, False,0,0) / (111.*1000.)
     dV_bardy   = center_finite_diff_n(v_bar, lat, False,0,0) / (111.*1000.)

    u2_dU_bardx = -1.*uu_pri_avg * dU_bardx
    uv_dU_bardy = -1.*uv_pri_avg * dU_bardy
    uv_dV_bardx = -1.*uv_pri_avg * dV_bardx
    v2_dV_bardy = -1.*vv_pri_avg * dV_bardy

    kp = u2_dU_bardx + uv_dU_bardy + uv_dV_bardx + v2_dV_bardy; + v2_v_bar

    copy_VarCoords(u_bar, kp) 
    return(kp)
end 
