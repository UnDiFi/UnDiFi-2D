      real(wp) uv,av,pv,rv
      real(wp) um,am,pm,rm
      real(wp) z1m,z2m,z3m,z4m
      real(wp) z1v,z2v,z3v,z4v
      real(wp) WsAvg,WsMax
      COMMON/SHOCKCOM/uv,av,pv,rv,um,am,pm,rm,                       &
     &z1v,z2v,z3v,z4v,z1m,z2m,z3m,z4m,WsAvg,WsMax

!    conditions upstream/downstream of the shock:
!    values are initialized in shockini
