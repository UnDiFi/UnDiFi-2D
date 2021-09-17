      ISLocalToGlobalMapping mapping(3)
      COMMON/COMMAP/mapping
C
C     mapping(1) is the mapping for the flow eqn. matrix
C     mapping(2) is the mapping for the turbulence eqn. matrix
C     mapping(3) is the mapping for the motion solver matrix
C
