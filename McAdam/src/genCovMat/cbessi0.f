C
C     compute complex modified Bessel functions I_0(z) of zeroth
C     order and first kind
C     KSM 28/5/2001
C
C     based on code by Shanjie Zhang and Jianming Jin
C     "Computation of Special Functions", 1996, John Wiley & Sons
C     routine CIK01, 
C     http://www.esrf.fr/computing/expg/libraries/smf/PROGRAMS/MCIK01.FOR
C
C     and Milton Abramovitz, Irene A. Stegun, 
C     "Handbook of Mathematical Functions", Dover Publications Inc.
C
C     PARAMS
C       Input :  C, D --- real and imaginary parts of complex arg
C       Output:  BC, BD --- real and imaginary parts of I_0(z)
C
      SUBROUTINE CDBESSI0(C, D, BC, BD)
      IMPLICIT NONE
      DOUBLE PRECISION C,D,BC,BD

      DOUBLE COMPLEX Z, CBI0
      DOUBLE COMPLEX CA, CR, Z2, Z1, ZR
      DOUBLE PRECISION PI, A0, A(12)
      INTEGER K, K0

      Z=DCMPLX(C,D)
      PI=3.141592653589793D0
      A0=CDABS(Z)
      IF (A0.EQ.0.0D0) THEN
         CBI0=(1.0D0,0.0D0)
         GOTO 15
      ENDIF

      Z2=Z*Z
      Z1=Z
C     make sure that | ARG Z | < 1/2 \pi; I_0 (z) = I_0 (-z)
      IF (DBLE(Z).LT.0.0) Z1=-Z

c     if MOD A < 18. use expansion
C     I_0 (Z) = 1 + (Z^2/4)/(1!)^2 + (Z^2/4)^2/(2!)^2 
C                 + (Z^2/4)^3/(3!)^2 + ...
      IF (A0.LE.18.0) THEN
         CBI0=(1.0D0,0.0D0)
         CR=(1.0D0,0.0D0)
         DO K=1,50
            CR=0.25D0*CR*Z2/(K*K)
            CBI0=CBI0+CR
            IF (CDABS(CR/CBI0).LT.1.0D-15) GOTO 15
         END DO
         
      ELSE
C     use asymptotic expansion for | ARG Z | < 1/2 \pi
C     I_o (z) ~ e^z / \sqrt{2 \pi z} { 1 + 1*1/(8z) + 3*3/(2!(8z)^2) 
C     + 1*1*3*3*5*5/(3!(8z)^3) + ...}
         DATA A/0.125D0,7.03125D-2,
     &        7.32421875D-2,1.1215209960938D-1,
     &        2.2710800170898D-1,5.7250142097473D-1,
     &        1.7277275025845D0,6.0740420012735D0,
     &        2.4380529699556D01,1.1001714026925D02,
     &        5.5133589612202D02,3.0380905109224D03/
         K0=12
         IF (A0.GE.35.0) K0=9
         IF (A0.GE.50.0) K0=7
         CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
         CBI0=(1.0D0,0.0D0)
         ZR=1.0D0/Z1
         DO K=1,K0
            CBI0=CBI0+A(K)*ZR**K
         END DO
         CBI0=CA*CBI0
      ENDIF
      
 15   BC = DBLE(CBI0)
      BD = DIMAG(CBI0)
      RETURN
      END

C     compute the modified complex Bessel function times the 
C     exponential of a given number
C     PARAMS
C       Input :  ARG -- real number by whose exponential the 
C                       Bessel function is multiplied
C       Input :  C, D --- real and imaginary parts of complex arg
C       Output:  BC, BD --- real and imaginary parts of I_0(z)
C
      SUBROUTINE CDEFACBESSI0(ARG, C, D, BC, BD)
      IMPLICIT NONE
      DOUBLE PRECISION ARG,C,D,BC,BD

      DOUBLE COMPLEX Z, CBI0
      DOUBLE COMPLEX CA, CR, Z2, Z1, ZR, ZA
      DOUBLE PRECISION PI, A0, A(12)
      INTEGER K, K0

      Z=DCMPLX(C,D)
      PI=4d0*atan(1d0)
      A0=CDABS(Z)
      IF (A0.EQ.0.0D0) THEN
         CBI0=(1.0D0,0.0D0)
         GOTO 16
      ENDIF

      Z2=Z*Z
      Z1=Z
C     make sure that | ARG Z | < 1/2 \pi; I_0 (z) = I_0 (-z)
      IF (DBLE(Z).LT.0.0) Z1=-Z

c     if MOD A < 18. use expansion 
C     I_0 (Z) = 1 + (Z^2/4)/(1!)^2 + (Z^2/4)^2/(2!)^2 
C                 + (Z^2/4)^3/(3!)^2 + ...
      IF (A0.LE.18.0) THEN
         CBI0=(1.0D0,0.0D0)
         CR=(1.0D0,0.0D0)
         DO K=1,50
            CR=0.25D0*CR*Z2/(K*K)
            CBI0=CBI0+CR
            IF (CDABS(CR/CBI0).LT.1.0D-15) GOTO 14
         END DO
 14      CBI0 = EXP(ARG) * CBI0
      ELSE
C     use asymptotic expansion for | ARG Z | < 1/2 \pi
C     I_o (z) ~ e^z / \sqrt{2 \pi z} { 1 + 1*1/(8z) + 3*3/(2!(8z)^2) 
C     + 1*1*3*3*5*5/(3!(8z)^3) + ...}
         DATA A/0.125D0,7.03125D-2,
     &        7.32421875D-2,1.1215209960938D-1,
     &        2.2710800170898D-1,5.7250142097473D-1,
     &        1.7277275025845D0,6.0740420012735D0,
     &        2.4380529699556D01,1.1001714026925D02,
     &        5.5133589612202D02,3.0380905109224D03/
         K0=12
         IF (A0.GE.35.0) K0=9
         IF (A0.GE.50.0) K0=7
         ZA=DCMPLX(ARG, 0.0D0)
         CA=CDEXP(ARG+Z1)/CDSQRT(2.0D0*PI*Z1)
         CBI0=(1.0D0,0.0D0)
         ZR=1.0D0/Z1
         DO K=1,K0
            CBI0=CBI0+A(K)*ZR**K
         END DO
         CBI0=CA*CBI0
      ENDIF
      
 16   BC = DBLE(CBI0)
      BD = DIMAG(CBI0)
      RETURN
      END

