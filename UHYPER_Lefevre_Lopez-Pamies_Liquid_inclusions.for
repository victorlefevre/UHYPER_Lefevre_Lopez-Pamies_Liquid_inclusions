!**********************************************************************
! Legal notice: UHYPER_Lefevre_Lopez-Pamies_Liquid_inclusions.for 
! (Windows)
!
! Copyright (C) 2020 Victor Lefèvre (victor.lefevre@northwestern.edu)
!                    Oscar Lopez-Pamies (pamies@illinois.edu)
!
! This ABAQUS UHYPER subroutine implements the hyperelastic energy  
! density derived in [1] for the macroscopic elastic response of
! isotropic and incompressible elastomers filled with liquid-like
! inclusions. The results applies to general non-percolative isotropic 
! distributions of incompressible inclusion whith vanishing shear
! stiffness. This result is valid for any choice of I1-based
! incompressible energy density characterizing the non-Gaussian 
! isotropic elastic response of the underlying elastomer. The present    
! subroutine is implemented for the choice of strain energy density
! proposed in [2].
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see https://www.gnu.org/licenses/
!
!**********************************************************************
! Usage:
!
! The subroutine is to be used as an incompressible USER hyperelastic 
! model with either 6 (arbitrary microstructures) or 5 (equiaxed
! particles) material properties, e.g.,
! *HYPERELASTIC, USER, TYPE=INCOMPRESSIBLE, PROPERTIES=6
! or
! *HYPERELASTIC, USER, TYPE=INCOMPRESSIBLE, PROPERTIES=5
! in the input (.inp) file.
!
! 1/ For arbitrary microstructures, the 6 materials properties required
! by the model to be input to the subroutine via the PROPS array are
! listed in the table below:
!
!  AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
!  ALPHA1  = PROPS(2)  ! EXPONENT #1 OF THE ELASTOMER
!  AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
!  ALPHA2  = PROPS(4)  ! EXPONENT #2 OF THE ELASTOMER 
!  AC      = PROPS(5)  ! VOLUME FRACTION OF LIQUID-LIKE INCLUSIONS
!  AMUT    = PROPS(6)  ! INITIAL SHEAR MODULUS OF THE 
!                      ! ELASTOMER COMPOSITE
!
! The two material parameters AMU1, AMU2 characterizing the elastic    
! behavior of the underlying elastomer are non-negative real numbers 
! (AMU1 >= 0, AMU2 >= 0) with strictly positive sum (AMU1 + AMU2 > 0). 
! The two exponents ALPHA1, ALPHA2 are non-zero real numbers 
! (ALPHA1 ≠ 0, ALPHA2 ≠ 0) leading to a strongly elliptic strain 
! energy (see eq. (22) in [2]). This is left to the user to check.
!
! The volume fractions of liquid-like inclusions (AC) must satisfy
! 0 <= AC <= 1.
!
! The initial shear modulus of the elastomer composite AMUT must be a 
! non-negative real number (AMUT >= 0).
!
! 2/ For microstructures comprising equiaxed liquid-like inclusions,  
! the 5 materials properties required by the model to be input to the
! subroutine via the PROPS array are listed in the table below:
!
!  AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
!  ALPHA1  = PROPS(2)  ! EXPONENT #1 OF THE ELASTOMER
!  AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
!  ALPHA2  = PROPS(4)  ! EXPONENT #2 OF THE ELASTOMER 
!  AC      = PROPS(5)  ! VOLUME FRACTION OF LIQUID-LIKE INCLUSIONS
!
! These 5 material properties are subjected to the same restrictions 
! listed above.
!
!**********************************************************************
! Additional information:
!
! This subroutine does not create solution-dependent state variables
! nor predefined field variables. 
!
! Please consult the ABAQUS Documentation for additional references 
! regarding the use of incompressible USER hyperelastic models with
! the UHYPER subroutine.
!
! Due the incompressible nature of this model, use of hybrid elements 
! is strongly recommended.
!
!**********************************************************************
! References:
!
! [1] Lefèvre, V., Lopez-Pamies, O. 2017. Nonlinear electroelastic 
!     deformations of dielectric elastomer composites: II — Non-Gaus-
!     sian elastic dielectrics. J. Mech. Phys. Solids  99, 438--470.
! [2] Lopez-Pamies, O., 2010. A new I1-based hyperelastic model for 
!     rubber elastic materials. C. R. Mec. 338, 3--11.
!
!**********************************************************************
!      
      SUBROUTINE UHYPER(BI1,BI2,AJ,U,UI1,UI2,UI3,TEMP,NOEL,
     1 CMNAME,INCMPFLAG,NUMSTATEV,STATEV,NUMFIELDV,FIELDV,
     2 FIELDVINC,NUMPROPS,PROPS)
!
      INCLUDE 'ABA_PARAM.INC'
#INCLUDE <SMAASPUSERSUBROUTINES.HDR>
!
      CHARACTER*80 CMNAME
      DIMENSION U(2),UI1(3),UI2(6),UI3(6),STATEV(*),FIELDV(*),
     1 FIELDVINC(*),PROPS(*)     
!
!     STDB_ABQERR AND GET_THREAD_ID INITIALIZATION
!
      DIMENSION INTV(1),REALV(4)
      CHARACTER*8 CHARV(1)
      CHARACTER*100 STRING1, STRING2, STRING3
      CHARACTER*300 STRING
!      
      INTEGER MYTHREADID  
!      
      INTV(1)=0
      REALV(1)=0.
      REALV(2)=0.
      REALV(3)=0.
      REALV(4)=0.
      CHARV(1)=''
!
      MYTHREADID = GET_THREAD_ID()
!
!     INPUT CHECKS
!
      IF (MYTHREADID.EQ.0) THEN
        IF (INCMPFLAG.EQ.0) THEN
          STRING1='INCOMPRESSIBILITY FLAG IS 0. THE MODEL IS INCOMPRES'
          STRING2='SIBLE. SET USER TYPE=INCOMPRESSIBLE.'
          STRING = TRIM(STRING1) // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF (NUMSTATEV.NE.0) THEN  
          INTV(1)=NUMSTATEV
          STRING1='RECEIVED REQUEST FOR %I SOLUTION-DEPENDENT STATE'
          STRING2=' VARIABLES. THE SUBROUTINE DOES NOT CREATE SOLUTION'
          STRING3='-DEPENDENT STATE VARIABLES.'
          STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF (NUMFIELDV.NE.0) THEN 
          INTV(1)=NUMFIELDV 
          STRING1='RECEIVED REQUEST FOR %I PREDEFINED FIELD  VARI'
          STRING2='ABLES. THE SUBROUTINE DOES NOT CREATE PREDEFINED'
          STRING3=' FIELD VARIABLES.'
          STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF ((NUMPROPS.NE.5).AND.(NUMPROPS.NE.6)) THEN    
          INTV(1)=NUMPROPS  
          STRING1='RECEIVED %I MATERIAL PROPERTIES. THE SUBROUTINE'
          STRING2=' REQUIRES EITHER 5 OR 6 MATERIAL PROPERTIES.'
          STRING = TRIM(STRING1) // TRIM(STRING2)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF
      END IF
!
!     MATERIAL PARAMETERS
!  
      AMU1    = PROPS(1)  ! PARAMETER #1 OF THE ELASTOMER
      ALPHA1  = PROPS(2)  ! I1 EXPONENT #1 OF THE ELASTOMER
      AMU2    = PROPS(3)  ! PARAMETER #2 OF THE ELASTOMER
      ALPHA2  = PROPS(4)  ! I1 EXPONENT #2 OF THE ELASTOMER 
      AC     = PROPS(5)  ! VOLUME FRACTION OF LIQUID-LIKE INCLUSIONS 
      IF (NUMPROPS.EQ.6) THEN
        AMUT  = PROPS(6)  ! INITIAL SHEAR MODULUS OF THE 
!                         ! ELASTOMER COMPOSITE
      END IF    
!
!     PARTIAL MATERIAL PARAMETERS CHECKS
!  
      IF (((AMU1.LT.0.).OR.(AMU2.LT.0.).OR.(AMU1+AMU2.LE.0.))
     1 .AND.(MYTHREADID.EQ.0)) THEN
        REALV(1)=AMU1      
        REALV(2)=AMU2      
        REALV(3)=AMU1+AMU2      
        STRING1='RECEIVED AMU1 = %R AND AMU2 = %R, AMU1 + AMU2 = %R.'
        STRING2=' THE PARAMETERS AMU1 AND AMU2 MUST BE NON-NEGATIVE'
        STRING3=' AND AMU1 + AMU2, MUST BE GREATER THAT ZERO.'
        STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      END IF           
!      
      IF (((ALPHA1.EQ.0.).OR.(ALPHA2.EQ.0.)).AND.
     1 (MYTHREADID.EQ.0)) THEN
        REALV(1)=ALPHA1      
        REALV(2)=ALPHA2      
        STRING1='RECEIVED ALPHA1 = %R AND ALPHA2 = %R.'
        STRING2=' THE EXPONENTS ALPHA1 AND ALPHA2 MUST BE NON-ZERO.'
        STRING = TRIM(STRING1) // TRIM(STRING2)
        CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
      END IF                     
!
      IF (AC.LT.0) THEN 
        IF (MYTHREADID.EQ.0) THEN
          REALV(1)=AC  
          STRING1='RECEIVED NEGATIVE AC = %R.' 
          STRING = TRIM(STRING1)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE
          CALL XIT
        END IF
      ELSE IF (AC.GE.1) THEN 
        IF (MYTHREADID.EQ.0) THEN
          REALV(1)=AC   
          STRING1='RECEIVED AC = %R GREATER OR EQUAL THAN 1.'
          STRING = TRIM(STRING1)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE
          CALL XIT
        END IF
      END IF    
!
      IF ((NUMPROPS.EQ.6).AND.(MYTHREADID.EQ.0)) THEN
        IF (AMUT.LT.0.) THEN
          REALV(1)=AMUT      
          STRING1='RECEIVED AMUT = %R.'
          STRING2=' THE INITIAL SHEAR MODULUS OF THE ELASTOMER COMPOSITE'
          STRING3=' MUST BE NON-NEGATIVE.'
          STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        ELSE IF (AMUT.GT.(AMU1+AMU2)) THEN  
          REALV(1)=AMUT       
          REALV(2)=AMU1+AMU2   
          STRING1='RECEIVED AMUT = %R. AND AMU1+AMU2 = %R.'
          STRING2=' THE INITIAL SHEAR MODULUS OF THE ELASTOMER COMPOSITE'
          STRING3=' MUST BE LOWER THAN THAT OF THE UNDERLYING ELASTOMER.'
          STRING = TRIM(STRING1) // TRIM(STRING2) // TRIM(STRING3)
          CALL STDB_ABQERR(-3,STRING,INTV,REALV,CHARV)
        END IF          
      END IF               
!      
      IF (NUMPROPS.EQ.6) THEN 
        ! ARBITRARY MICROSTRUCTURE [1]
        ALC = AMUT/((AMU1+AMU2)*(1.-AC))
      ELSEIF (NUMPROPS.EQ.5) THEN
        ! EQUIAXED LIQUID-LIKE INCLUSIONS [1]
        ALC = (1.-AC)**(2.0/3.0)
      END IF
!      
!     ‘‘AMPLIFIED’’ STRAIN MEASURE 
!      
      AII1 = ALC*(BI1-3.)+3.
!      
!     NON-GAUSSIAN HYPERELASTIC MODEL FOR THE ELASTOMER [2] 
!     AND DERIVATIVES
!      
      AP1 = AMU1*3.**(1.-ALPHA1)*0.5
      AP2 = AMU2*3.**(1.-ALPHA2)*0.5   
!      
      PSI = AP1/ALPHA1*(AII1**ALPHA1-3.**ALPHA1)+
     1      AP2/ALPHA2*(AII1**ALPHA2-3.**ALPHA2)
!     
      DPSI = AP1*AII1**(ALPHA1-1.)+AP2*AII1**(ALPHA2-1.)
!     
      DDPSI = AP1*(ALPHA1-1.)*AII1**(ALPHA1-2.)+
     1        AP2*(ALPHA2-1.)*AII1**(ALPHA2-2.)  
!
!     MACROSCOPIC STRAIN ENERGY DENSITY FUNCTION FOR
!     FILLED ELASTOMER [1,3,4]
!      
      U(1) = (1.-AC)*PSI
      U(2) = 0.
!
!     FIRST PARTIAL DERIVATIVES
!
      UI1(1) = (1.-AC)*ALC*DPSI
      UI1(2) = 0.
      UI1(3) = 0.
!
!     SECOND PARTIAL DERIVATIVES
!
      UI2(1) = (1.-AC)*ALC**2.*DDPSI
      UI2(2) = 0.
      UI2(3) = 0.
      UI2(4) = 0.
      UI2(5) = 0.
      UI2(6) = 0.
!
!     THIRD PARTIAL DERIVATIVES
!
      UI3(1) = 0.
      UI3(2) = 0.
      UI3(3) = 0.
      UI3(4) = 0.
      UI3(5) = 0.
      UI3(6) = 0.   
!      
      RETURN
      END
!
!**********************************************************************