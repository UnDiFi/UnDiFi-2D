      LOGICAL       TIMEIMPL,PICARD,NEWTON
      INTEGER       IGLOB,CFL_SER,CFL_EXP
      DOUBLE PRECISION PE
      COMMON/IMPSOL/TIMEIMPL,PICARD,NEWTON,IGLOB
      PARAMETER(CFL_SER=1,CFL_EXP=2,PE=30.d0)
C
C     TIMEIMPL = .TRUE.  ---->     use implicit time-stepping
C     PICARD = .TRUE.  ---->     use Picard linearization
C     NEWTON = .TRUE.  ---->     use F.D. Newton linearization
C     IGLOB  = 0 (.TRUE.)  ---> use global timestepping
C
C     [default]:
C
C     TIMEIMPL = .TRUE.
C     PICARD = .TRUE. 
C     NEWTON = .FALSE. 
C     IGLOB  = 1 (.FALSE.)  ---> use Local timestepping
C     option -explicit sets TIMEIMPL = .FALSE.
C     option -Newton   sets PICARD = .FALSE. .AND. NEWTON = .TRUE.
