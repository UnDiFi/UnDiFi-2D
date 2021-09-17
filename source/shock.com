      DOUBLE PRECISION uv,av,pv,rv
      DOUBLE PRECISION um,am,pm,rm
      DOUBLE PRECISION z1m,z2m,z3m,z4m
      DOUBLE PRECISION z1v,z2v,z3v,z4v
      DOUBLE PRECISION WsAvg,WsMax
      COMMON/SHOCKCOM/uv,av,pv,rv,um,am,pm,rm,
     &z1v,z2v,z3v,z4v,z1m,z2m,z3m,z4m,WsAvg,WsMax

c    conditions upstream/downstream of the shock:
c    values are initialized in shockini
