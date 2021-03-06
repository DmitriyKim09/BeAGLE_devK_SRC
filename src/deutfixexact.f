      SUBROUTINE DEUTFIXEXACT(NU,Q2,MDEUT,SINDEX)
C
C     2020-03-09 added by Kong Tu
C
      IMPLICIT NONE
      DOUBLE PRECISION NU, Q2, MDEUT

      include 'beagle.inc'
C      include "py6strf.inc"   ! Temporary! Just use for debug output

C      include 'pythia.inc' - conflicts with IMPLICIT NONE
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.0D0)

CCc...added by liang & Mark to include pythia energy loss datas
C      double precision PAUX, DPF
C      COMMON /PFAUX/ PAUX(4), DPF(4)

* event flag
      COMMON /DTEVNO/ NEVENT,ICASCA 
      integer NEVENT, ICASCA

C Local
      DOUBLE PRECISION MJ, MP, MN
      INTEGER NDIM,MAXPRTS
      PARAMETER (NDIM=5)
      PARAMETER (MAXPRTS=20)
      PARAMETER (MJ=3.09688D0)
      PARAMETER (MP=0.93827D0)
      PARAMETER (MN=0.93957D0)
      DOUBLE PRECISION PPL(MAXPRTS,5),PSPEC(NDIM)
      DOUBLE PRECISION QZKZ,NUMN
      DOUBLE PRECISION PX,PY,PZ
      DOUBLE PRECISION JX,JY,JZ
      DOUBLE PRECISION MSTRU,MSPEC

      INTEGER JNDEX,PNDEX
      INTEGER NPRTNS,NLSCAT,IDIM,ITRK,JTRK
      INTEGER INDXP(MAXPRTS)
      INTEGER SINDEX

      IF (IOULEV(4).GE.2 .AND. NEVENT.LE.IOULEV(5)) THEN
         WRITE (*,*) 'DEUTFIXEXACT: About to fix e+D kinematics exact'
         CALL PYLIST(2)
         WRITE (*,*)
      ENDIF

      NLSCAT = 0
      NPRTNS = 0
      DO IDIM=1,NDIM
         PSPEC(IDIM)=ZERO     
      ENDDO
      DO ITRK=1,N
         IF(K(ITRK,1).EQ.1 .OR. K(ITRK,1).EQ.2) THEN
            IF ( (ABS(K(ITRK,2)).EQ.11 .OR. ABS(K(ITRK,2)).EQ.13) .AND.
     &           K(ITRK,3).EQ.3) THEN
               NLSCAT = NLSCAT+1
            ELSE
               IF( ITRK .NE. SINDEX ) THEN
                  NPRTNS=NPRTNS+1
                  IF (NPRTNS.GT.MAXPRTS) 
     &              STOP('DEUTFIXEXACT: FATAL ERROR. Too many partons')
                  INDXP(NPRTNS)=ITRK
                  DO IDIM=1,NDIM
                    PPL(NPRTNS,IDIM)=P(ITRK,IDIM)
                  ENDDO
C     identify spectator     
               ELSE
                  DO IDIM=1,NDIM
                    PSPEC(IDIM)=P(ITRK,IDIM)
                  ENDDO  
                  MSPEC = P(ITRK,5)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF (NLSCAT.NE.1) 
     &     STOP "ERROR! BAD EVENT CONFIG. Scattered leptons .ne. 1"
      IF (NPRTNS.LT.2)
     &     STOP "ERROR! BAD EVENT CONFIG. Fewer than two particles"

C     a quick and dirty way of getting struck mass
      IF( MSPEC < 0.939D0 ) THEN
        MSTRU = MN
      ELSE
        MSTRU = MP
      ENDIF

      QZKZ = SQRT(NU**2+Q2) - PSPEC(3)
      NUMN = NU - PSPEC(4)

C     Initialize px,py,jx,jy
      PX = 0.0D0
      PY = 0.0D0
      JX = 0.0D0
      JY = 0.0D0

C     Find Jpsi and struck nucleon
      DO JTRK=1,NPRTNS
        IF( ABS(PPL(JTRK,5)-MJ) .LT. 0.0001D0 ) THEN
          JX = PPL(JTRK,1)
          JY = PPL(JTRK,2)
          JNDEX = JTRK
        ENDIF
        IF( ABS(PPL(JTRK,5)-MSTRU) .LT. 0.0001D0 ) THEN
          PX = PPL(JTRK,1)
          PY = PPL(JTRK,2)
          PNDEX = JTRK
        ENDIF
      ENDDO

C     solution for struck nucleon z momentum

      PZ = (QZKZ*(-JX**2 - JY**2 - MJ**2 + MSTRU**2 + (MDEUT + NUMN)**2 
     &   + PX**2 + PY**2 - 
     & QZKZ**2) - Sqrt((MDEUT + NUMN)**2*
     & (JX**4 + JY**4 + MDEUT**4 - 2*MDEUT**2*MJ**2 + MJ**4 - 
     &  2*MDEUT**2*MSTRU**2 - 
     & 2*MJ**2*MSTRU**2 + MSTRU**4 + 4*MDEUT**3*NUMN - 
     &  4*MDEUT*MJ**2*NUMN - 
     & 4*MDEUT*MSTRU**2*NUMN + 6*MDEUT**2*NUMN**2 - 2*MJ**2*NUMN**2 - 
     & 2*MSTRU**2*NUMN**2 + 4*MDEUT*NUMN**3 + NUMN**4 - 
     &  2*MDEUT**2*PX**2 - 
     & 2*MJ**2*PX**2 + 2*MSTRU**2*PX**2 - 4*MDEUT*NUMN*PX**2 - 
     &  2*NUMN**2*PX**2 + 
     & PX**4 - 2*MDEUT**2*PY**2 - 2*MJ**2*PY**2 + 2*MSTRU**2*PY**2 - 
     & 4*MDEUT*NUMN*PY**2 - 2*NUMN**2*PY**2 + 2*PX**2*PY**2 + PY**4 + 
     & 2*(MJ**2 + MSTRU**2 - (MDEUT + NUMN)**2 + PX**2 + PY**2)*QZKZ**2 
     &  + QZKZ**4 - 
     & 2*JY**2*(-MJ**2 + MSTRU**2 + (MDEUT + NUMN)**2 + PX**2 + PY**2 
     &  - QZKZ**2) + 
     & 2*JX**2*(JY**2 + MJ**2 - MSTRU**2 - (MDEUT + NUMN)**2 - 
     & PX**2 - PY**2 + 
     &    QZKZ**2))))/(2.*(MDEUT + NUMN - QZKZ)*(MDEUT + NUMN + QZKZ))

C     solution for Jpsi z momentum

      JZ = (QZKZ*(JX**2 + JY**2 + MJ**2 - MSTRU**2 + (MDEUT + NUMN)**2  
     & - PX**2 - PY**2 - QZKZ**2) + Sqrt((MDEUT + NUMN)**2*
     & (JX**4 + JY**4 + MDEUT**4 - 2*MDEUT**2*MJ**2 
     & + MJ**4 - 2*MDEUT**2*MSTRU**2 - 
     & 2*MJ**2*MSTRU**2 + MSTRU**4 + 4*MDEUT**3*NUMN 
     & - 4*MDEUT*MJ**2*NUMN - 
     & 4*MDEUT*MSTRU**2*NUMN + 6*MDEUT**2*NUMN**2 - 2*MJ**2*NUMN**2 - 
     & 2*MSTRU**2*NUMN**2 + 4*MDEUT*NUMN**3 + NUMN**4 - 
     & 2*MDEUT**2*PX**2 - 
     & 2*MJ**2*PX**2 + 2*MSTRU**2*PX**2 - 4*MDEUT*NUMN*PX**2 
     & - 2*NUMN**2*PX**2 + 
     & PX**4 - 2*MDEUT**2*PY**2 - 2*MJ**2*PY**2 + 2*MSTRU**2*PY**2 - 
     & 4*MDEUT*NUMN*PY**2 - 2*NUMN**2*PY**2 + 2*PX**2*PY**2 + PY**4 + 
     & 2*(MJ**2 + MSTRU**2 - (MDEUT + NUMN)**2 + PX**2 + PY**2)*QZKZ**2 
     & + QZKZ**4 - 2*JY**2*(-MJ**2 + MSTRU**2
     &  + (MDEUT + NUMN)**2 + PX**2 + PY**2 - QZKZ**2) + 
     & 2*JX**2*(JY**2 + MJ**2 - MSTRU**2 - (MDEUT + NUMN)**2 - PX**2  
     & - PY**2 + QZKZ**2))))/(2.*(MDEUT + NUMN - QZKZ)*
     & (MDEUT + NUMN + QZKZ))

      PPL(JNDEX,3) = JZ
      PPL(JNDEX,4) = SQRT(JX**2+JY**2+JZ**2+P(INDXP(JNDEX),5)**2)
      PPL(PNDEX,3) = PZ
      PPL(PNDEX,4) = SQRT(PX**2+PY**2+PZ**2+P(INDXP(PNDEX),5)**2)

C     only modify Jpsi and struck nucleon

      DO IDIM=1,NDIM
        P(INDXP(JNDEX),IDIM)=PPL(JNDEX,IDIM)
        P(INDXP(PNDEX),IDIM)=PPL(PNDEX,IDIM)
      ENDDO
      
      RETURN
      END