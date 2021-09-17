      DOUBLE PRECISION DELT,GAMT,GAMTBAK,THETAT,DTVOL,TCOEF(-1:1),
     2                 ALFALE,XYZ_C(3),OPITCH(3),APITCH(3)
      INTEGER ITSTEP,NTIMLEVS,MMTYPE,ITIM,IALE
      LOGICAL LTIME,DUALTS,LFOOLD,LALE,ALE_MOVE_GRID,ALE_READ_GRID,
     2        ALE_LAPLACE_GRID,ALE_READ_GRIDVEL,CHAR_TIMESTEPPING,GCLCHK
      COMMON/R8TIME/DELT,GAMT,GAMTBAK,THETAT,DTVOL,TCOEF,ALFALE,XYZ_C,
     2              OPITCH,APITCH
      COMMON/I4TIME/ITSTEP,NTIMLEVS,MMTYPE,ITIM,IALE
      COMMON/L4TIME/LTIME,DUALTS,LFOOLD,LALE,ALE_MOVE_GRID,
     2              ALE_READ_GRID,ALE_LAPLACE_GRID,ALE_READ_GRIDVEL,
     3              CHAR_TIMESTEPPING,GCLCHK
C
!>     DELT is the time step size for time accurate calculations (-time_step_size)
!>     GAMT weights the contribution of the three time levels: (n+1), n, (n-1)
!>     GAMTBAK is a backup copy of GAMT
!>     THETAT weights the contribution of the spatial residuals at time level (n+1) and n
!>     DTVOL 
!>     \todo DTVOL is a local parameter, should be rather placed elsewhere
!>     TCOEF(-1:1) coefficients of the 3 time levels scheme
!>     ALFALE weights used to determine the mesh configurations that guarantees GCL
!>     XYZ_C coordinates of the point about which rotation/deformation is performed when moving the grid
!>     OPITCH dimensionless angular velocity of the pitching motion
!>     APITCH pitching amplitude velocity of the pitching motion
!>
!>     ITSTEP is the number of time steps (-nof_time_steps)
!>     NTIMLEVS is the actual nof of time levels; can be two: (n+1) and n or three: (n+1), n, (n-1)
!>     MMTYPE is the scheme used for the Mass Matrix, see time.h 
!>     ITIM is the current time step, i.e. we are moving from n to n+1
!>     IALE control the function used to move the grid (or its boundaries)
!>     LTIME is set to .TRUE. for time accurate calculations; activated using:
!>
!>     DUALTS is set to .TRUE. when dual time stepping is chosen (-dual_ts [Y/N])
!>     LFOOLD should be set to .TRUE. when reading time level n-1 from a datafile
!>     LALE     = .TRUE. when doing arbitrary eulerian lagrangian 
!>     \verbatim
!>     -ale
!>     \endverbatim
!>     ALE_MOVE_GRID when .TRUE. the new mesh is only moved
!>     ALE_LAPLACE_GRID when .TRUE. the new mesh is obtained by solving Laplace's equation for the grid velocities
!>     ALE_READ_GRID when .TRUE. the new mesh is read from file (currently un-used)
!>     ALE_READ_GRIDVEL when .TRUE. the new mesh is read from file
!>     CHAR_TIMESTEPPING is set to .TRUE. when characteristic time stepping is chosen
!>     GCLCHK when set to .TRUE. checks whether the GCL is verified on a per cell basis; activated using
!>     \verbatim
!>     -ale_check_gcl
!>     \endverbatim
!>
