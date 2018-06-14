      DOUBLE PRECISION FUNCTION ALPHAS_ellis(Q,AMZ,NLOOP)
c     Evaluation of strong coupling constant alpha_S
c     Author: R.K. Ellis

c     q -- scale at which alpha_s is to be evaluated
c     amz -- value of alpha_s at the mass of the Z-boson
c     nloop -- the number of loops (1,2, or 3) at which beta
c     function is evaluated to determine running.
c     the values of the cmass and the bmass should be set
c     in common block qmass.

      IMPLICIT NONE
      DOUBLE PRECISION Q,T,AMZ,AMZ0,AMB,AMC,ZMASS,BMASS,CMASS,AS_OUT
      INTEGER NLOOP,NLOOP0,NF3,NF4,NF5
      PARAMETER(ZMASS=91.188D0)
      PARAMETER(NF5=5,NF4=4,NF3=3)
      COMMON/QMASS/CMASS,BMASS
      SAVE AMZ0,NLOOP0,AMB,AMC
      DATA AMZ0,NLOOP0/0D0,0/

      IF (Q .LE. 0D0) THEN
         WRITE(6,*) 'q .le. 0 in alphas'
         WRITE(6,*) 'q= ',Q
         STOP
      ENDIF
      IF (AMZ .LE. 0D0) THEN
         WRITE(6,*) 'amz .le. 0 in alphas',AMZ
         STOP
      ENDIF
      IF (CMASS .LE. 0.3D0) THEN
         WRITE(6,*) 'cmass .le. 0.3GeV in alphas',CMASS
         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
         WRITE(6,*) 'continue with cmass=1.5GeV'
         CMASS=1.5D0
      ENDIF
      IF (BMASS .LE. 0D0) THEN
         WRITE(6,*) 'bmass .le. 0 in alphas',BMASS
         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
         WRITE(6,*) 'continue with bmass=5.0GeV'
         BMASS=5D0
      ENDIF
c--- establish value of coupling at b- and c-mass and save
      IF ((AMZ .NE. AMZ0) .OR. (NLOOP .NE. NLOOP0)) THEN
         AMZ0=AMZ
         NLOOP0=NLOOP
         T=2D0*DLOG(BMASS/ZMASS)
         CALL NEWTON1(T,AMZ,AMB,NLOOP,NF5)
         T=2D0*DLOG(CMASS/BMASS)
         CALL NEWTON1(T,AMB,AMC,NLOOP,NF4)
      ENDIF

c--- evaluate strong coupling at scale q
      IF (Q  .LT. BMASS) THEN
           IF (Q  .LT. CMASS) THEN
             T=2D0*DLOG(Q/CMASS)
             CALL NEWTON1(T,AMC,AS_OUT,NLOOP,NF3)
           ELSE
             T=2D0*DLOG(Q/BMASS)
             CALL NEWTON1(T,AMB,AS_OUT,NLOOP,NF4)
           ENDIF
      ELSE
      T=2D0*DLOG(Q/ZMASS)
      CALL NEWTON1(T,AMZ,AS_OUT,NLOOP,NF5)
      ENDIF
      ALPHAS_ellis=AS_OUT
      RETURN
      END

      SUBROUTINE SETQMASS(CMASS_IN, BMASS_IN)
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: CMASS_IN, BMASS_IN
      COMMON/QMASS/CMASS,BMASS
      DOUBLE PRECISION CMASS, BMASS
C
      CMASS = CMASS_IN
      BMASS = BMASS_IN
C
      END SUBROUTINE SETQMASS
