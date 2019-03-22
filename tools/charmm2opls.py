#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = range

import	re
import	math
import	collections


if( len( sys.argv ) != 3 ):
	print( "%s Topology(RTF) Parameters(PRM)"%( sys.argv[0] ) )
	sys.exit( 1 )


opls_type = """MM_Definitions OPLS_AA 1.0

Types
H       1   0.00000   0.00000
HA      1   2.42000   0.03000
HC      1   2.50000   0.03000
HO      1   0.00000   0.00000
HS      1   0.00000   0.00000
H3      1   0.00000   0.00000
C       6   3.75000   0.10500
C*      6   3.55000   0.07000
CA      6   3.55000   0.07000
CB      6   3.55000   0.07000
CCA     6   3.75000   0.10500
CN      6   3.55000   0.07000
CAA     6   2.25000   0.05000
CP      6   3.55000   0.07000
CT      6   3.50000   0.06600
CTP     6   3.50000   0.06600
CV      6   3.55000   0.07000
CW      6   3.55000   0.07000
N       7   3.25000   0.17000
NA      7   3.25000   0.17000
NB      7   3.25000   0.17000
N2      7   3.25000   0.17000
N3      7   3.25000   0.17000
O       8   2.96000   0.21000
OH      8   3.12000   0.17000
OHT     8   3.07000   0.17000
OHC     8   3.00000   0.17000
O2      8   2.96000   0.21000
OS      8   3.00000   0.17000
S      16   3.55000   0.25000
SH     16   3.55000   0.25000
OW      8   3.15061   0.15210 
HW      1   0.00000   0.00000
CL     17   4.41724   0.11779
BR     35   4.62376   0.09000
SOD    11   1.89744   1.60714
ZN     30   1.94215   0.25000
! ---------------------------
"""

opls_resi = """
Electrostatics Scale 0.5
Lennard_Jones  Scale 0.5
Units kcal/mole

Residues

Residue TIP3
  3  3  0
OH2      OW      -0.834
H1       HW       0.417
H2       HW       0.417

OH2 H1 ; OH2 H2 ; H1 H2

Residue HOH
  3  3  0
OH2      OW      -0.834
H1       HW       0.417
H2       HW       0.417

OH2 H1 ; OH2 H2 ; H1 H2

Residue ALA
 10 11  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.18
HB1   HC     0.06
HB2   HC     0.06
HB3   HC     0.06
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA  ; CA   HA ; CA   CB ; CA   C
CB   HB1 ; CB   HB2 ; CB   HB3 ; C    O  ; C    +R

-R   CA   N    H ; CA   +R   C    O                                               

Residue ARG
 24 25 6
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    CT    -0.05
HG1   HC     0.06
HG2   HC     0.06
CD    CT     0.19
HD1   HC     0.06
HD2   HC     0.06
NE    N2    -0.70
HE    H3     0.44
CZ    CAA    0.64
NH1   N2    -0.80
HH11  H3     0.46
HH12  H3     0.46
NH2   N2    -0.80
HH21  H3     0.46
HH22  H3     0.46
C     C      0.50
O     O     -0.50

-R   N   ; N    H    ; N    CA   ; CA   HA   ; CA   CB   ; CA   C
CB   HB1 ; CB   HB2  ; CB   CG   ; CG   HG1  ; CG   HG2  ; CG   CD
CD   HD1 ; CD   HD2  ; CD   NE   ; NE   HE   ; NE   CZ   ; CZ   NH1
CZ   NH2 ; NH1  HH11 ; NH1  HH12 ; NH2  HH21 ; NH2  HH22 ; C    O
C    +R

-R   CA   N    H    ; CA   +R   C    O                                               
NE   NH1  CZ   NH2  ; CZ   CD   NE   HE                                              
CZ   HH11 NH1  HH12 ; CZ   HH21 NH2  HH22

Residue ASN
 14 15  4
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    C      0.50
OD1   O     -0.50
ND2   N     -0.76
HD21  H      0.38
HD22  H      0.38
C     C      0.50
O     O     -0.50

-R   N    ; N    H   ; N    CA ; CA   HA  ; CA   CB  ; CA   C
CB   HB1  ; CB   HB2 ; CB   CG ; CG   OD1 ; CG   ND2 ; ND2  HD21
ND2  HD22 ; C    O   ; C    +R

-R   CA   N    H   ; CA   +R   C    O                                               
CB   ND2  CG   OD1 ; CG   HD21 ND2  HD22

Residue ASP
 12 13  3
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.22
HB1   HC     0.06
HB2   HC     0.06
CG    CCA    0.70
OD1   O2    -0.80
OD2   O2    -0.80
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG ; CG   OD1 ; CG   OD2 ; C    O
C    +R

-R   CA   N    H ; CA   +R   C    O ; CB   OD1  CG   OD2 

! -- Residue ASH

Residue CYS
 11 12  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT     0.06
HB1   HC     0.06
HB2   HC     0.06
SG    SH    -0.435
HG    HS     0.255
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA ; CA   HA ; CA   CB ; CA   C
CB   HB1 ; CB   HB2 ; CB   SG ; SG   HG ; C    O  ; C    +R

-R   CA   N    H ; CA   +R   C    O

Residue GLN
 17 18  4
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    CT    -0.12
HG1   HC     0.06
HG2   HC     0.06
CD    C      0.50
OE1   O     -0.50
NE2   N     -0.76
HE21  H      0.38
HE22  H      0.38
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA   ; CA   HA   ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG   ; CG   HG1  ; CG   HG2 ; CG   CD
CD   OE1 ; CD   NE2 ; NE2  HE21 ; NE2  HE22 ; C    O   ; C    +R

-R   CA   N    H   ; CA   +R   C    O                                               
CG   NE2  CD   OE1 ; CD   HE21 NE2  HE22

Residue GLU
 15 16  3
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    CT    -0.22
HG1   HC     0.06
HG2   HC     0.06
CD    CCA    0.70
OE1   O2    -0.80
OE2   O2    -0.80
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG ; CG   HG1 ; CG   HG2 ; CG   CD
CD   OE1 ; CD   OE2 ; C    O  ; C    +R

-R   CA   N    H ; CA   +R   C    O ; CG   OE1  CD   OE2                                             

! -- Residue GLH

Residue GLY                                                        
  7  8  2
N     N     -0.50
H     H      0.30
CA    CT     0.08
HA1   HC     0.06
HA2   HC     0.06
C     C      0.50
O     O     -0.50

-R   N   ; N    H ; N    CA ; CA   HA1
CA   HA2 ; CA   C ; C    O  ; C    +R

-R   CA   N    H ; CA   +R   C    O                                               

Residue HSD
 17 19  6
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.005
HB1   HC     0.06
HB2   HC     0.06
CG    CW     0.015
ND1   NA    -0.570
HD1   H      0.420
CE1   CP     0.295
HE1   HA     0.115
NE2   NB    -0.490
CD2   CV    -0.015
HD2   HA     0.115
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA  ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG  ; CG   ND1 ; ND1  HD1 ; ND1  CE1
CE1  HE1 ; CE1  NE2 ; NE2  CD2 ; CD2  HD2 ; CD2  CG  ; C    O
C    +R

-R   CA   N    H   ; CA   +R   C    O   ; CE1  CG   ND1  HD1                                             
CG   HD2  CD2  NE2 ; HE1  ND1  CE1  NE2 ; CB   CD2  CG   ND1

Residue HSE
 17 19  6
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.005
HB1   HC     0.06
HB2   HC     0.06
CG    CV    -0.015
ND1   NB    -0.490
CE1   CP     0.295
HE1   HA     0.115
NE2   NA    -0.570
HE2   H      0.420
CD2   CW     0.015
HD2   HA     0.115
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA  ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG  ; CG   ND1 ; ND1  CE1 ;CE1  HE1
CE1  NE2 ; NE2  HE2 ; NE2  CD2 ; CD2  HD2 ; CD2  CG  ; C    O
C    +R

-R   CA   N    H   ; CA   +R   C    O   ; CE1  CD2  NE2  HE2                                             
CG   HD2  CD2  NE2 ; HE1  NE2  CE1  ND1 ; CB   CD2  CG   ND1

Residue HSP
 18 20  7
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.005
HB1   HC     0.06
HB2   HC     0.06
CG    CW     0.215
ND1   NA    -0.540
HD1   H      0.460
CE1   CP     0.385
HE1   HA     0.115
NE2   NA    -0.540
HE2   H      0.460
CD2   CW     0.215
HD2   HA     0.115
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA  ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG  ; CG   ND1 ; ND1  HD1 ; ND1  CE1
CE1  HE1 ; CE1  NE2 ; NE2  HE2 ; NE2  CD2 ; CD2  HD2 ; CD2  CG  
C    O   ; C    +R

-R   CA   N    H   ; CA   +R   C    O   ; CE1  CG   ND1  HD1                                             
CG   HD2  CD2  NE2 ; HE1  ND1  CE1  NE2 ; CE1  CD2  NE2  HE2
CB   CD2  CG   ND1

Residue ILE                                                 
 19 20  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.06
HB    HC     0.06
CG2   CT    -0.18
HG21  HC     0.06
HG22  HC     0.06
HG23  HC     0.06
CG1   CT    -0.12
HG11  HC     0.06
HG12  HC     0.06
CD    CT    -0.18
HD1   HC     0.06
HD2   HC     0.06
HD3   HC     0.06
C     C      0.50
O     O     -0.50

-R   N    ; N    H    ; N    CA   ; CA   HA   ; CA   CB   ; CA   C
CB   HB   ; CB   CG1  ; CB   CG2  ; CG1  HG11 ; CG1  HG12 ; CG1  CD 
CG2  HG21 ; CG2  HG22 ; CG2  HG23 ; CD   HD1  ; CD   HD2  ; CD   HD3
C    O    ; C    +R

-R   CA   N    H ; CA   +R   C    O                                               

Residue LEU
 19 20  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    CT    -0.06
HG    HC     0.06
CD1   CT    -0.18
HD11  HC     0.06
HD12  HC     0.06
HD13  HC     0.06
CD2   CT    -0.18
HD21  HC     0.06
HD22  HC     0.06
HD23  HC     0.06
C     C      0.50
O     O     -0.50

-R   N    ; N    H    ; N    CA   ; CA   HA   ; CA   CB   ; CA   C
CB   HB1  ; CB   HB2  ; CB   CG   ; CG   HG   ; CG   CD1  ; CG   CD2
CD1  HD11 ; CD1  HD12 ; CD1  HD13 ; CD2  HD21 ; CD2  HD22 ; CD2  HD23
C    O    ; C    +R

-R   CA   N    H ; CA   +R   C    O

Residue LYS
 22 23  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    CT    -0.12
HG1   HC     0.06
HG2   HC     0.06
CD    CT    -0.12
HD1   HC     0.06
HD2   HC     0.06
CE    CT     0.19
HE1   HC     0.06
HE2   HC     0.06
NZ    N3    -0.30
HZ1   H3     0.33
HZ2   H3     0.33
HZ3   H3     0.33
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA  ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG  ; CG   HG1 ; CG   HG2 ; CG   CD
CD   HD1 ; CD   HD2 ; CD   CE  ; CE   HE1 ; CE   HE2 ; CE   NZ
NZ   HZ1 ; NZ   HZ2 ; NZ   HZ3 ; C    O   ; C    +R

-R   CA   N    H ; CA   +R   C    O

! -- Residue LYN

Residue MET
 17 18  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    CT     0.0975
HG1   HC     0.06
HG2   HC     0.06
SD    S     -0.435
CE    CT     0.0375
HE1   HC     0.06
HE2   HC     0.06
HE3   HC     0.06
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA  ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG  ; CG   HG1 ; CG   HG2 ; CG   SD
SD   CE  ; CE   HE1 ; CE   HE2 ; CE   HE3 ; C    O   ; C    +R

-R   CA   N    H ; CA   +R   C    O

Residue PHE
 20 22  8
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.005
HB1   HC     0.06
HB2   HC     0.06
CG    CA    -0.115
CD1   CA    -0.115
HD1   HA     0.115
CE1   CA    -0.115
HE1   HA     0.115
CZ    CA    -0.115
HZ    HA     0.115
CE2   CA    -0.115
HE2   HA     0.115
CD2   CA    -0.115
HD2   HA     0.115
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG ; CG   CD1 ; CD1  HD1 ; CD1  CE1
CE1  HE1 ; CE1  CZ  ; CZ   HZ ; CZ   CE2 ; CE2  HE2 ; CE2  CD2
CD2  HD2 ; CD2  CG  ; C    O  ; C    +R

-R   CA   N    H  ; CA   +R   C    O   ; CG   CE2  CD2  HD2 ; CD2  CZ   CE2  HE2
CE1  CE2  CZ   HZ ; CD1  CZ   CE1  HE1 ; CG   CE1  CD1  HD1 ; CD1  CD2  CG   CB

Residue PRO
 14 16  2
N     N     -0.14
CD    CTP   -0.05
HD1   HC     0.06
HD2   HC     0.06
CG    CT    -0.12
HG1   HC     0.06
HG2   HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CA    CT     0.01
HA    HC     0.06
C     C      0.50
O     O     -0.50

-R   N   ; N    CD  ; N    CA ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG ; CG   HG1 ; CG   HG2 ; CG   CD
CD   HD1 ; CD   HD2 ; C    O  ; C    +R

CA   +R   C    O ; -R   CA   N    CD

Residue SER
 11 12  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT     0.145
HB1   HC     0.06
HB2   HC     0.06
OG    OH    -0.683
HG    HO     0.418
C     C      0.50
O     O     -0.50

-R   N   ; N    H   ; N    CA ; CA   HA ; CA   CB ; CA   C
CB   HB1 ; CB   HB2 ; CB   OG ; OG   HG ; C    O  ; C    +R

-R   CA   N    H ; CA   +R   C    O                                               

Residue THR
 14 15  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT     0.205
HB    HC     0.06
CG2   CT    -0.18
HG21  HC     0.06
HG22  HC     0.06
HG23  HC     0.06
OG1   OH    -0.683
HG1   HO     0.418
C     C      0.50
O     O     -0.50

-R   N    ; N    H   ; N    CA  ; CA   HA  ; CA   CB   ; CA   C
CB   HB   ; CB   OG1 ; CB   CG2 ; OG1  HG1 ; CG2  HG21 ; CG2  HG22
CG2  HG23 ; C    O   ; C    +R

-R   CA   N    H ; CA   +R   C    O                                               

Residue TRP
 24 27  9
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.12
HB1   HC     0.06
HB2   HC     0.06
CG    C*     0.075
CD1   CA    -0.115
HD1   HA     0.115
NE1   NA    -0.570
HE1   H      0.420
CE2   CN     0.130
CZ2   CA    -0.115
HZ2   HA     0.115
CH2   CA    -0.115
HH2   HA     0.115
CZ3   CA    -0.115
HZ3   HA     0.115
CE3   CA    -0.115
HE3   HA     0.115
CD2   CB    -0.055
C     C      0.50
O     O     -0.50
                                                                
-R   N   ; N    H   ; N    CA  ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG  ; CG   CD1 ; CG   CD2 ; CD1  HD1
CD1  NE1 ; NE1  HE1 ; NE1  CE2 ; CD2  CE2 ; CD2  CE3 ; CE3  HE3
CE2  CZ2 ; CZ2  HZ2 ; CZ2  CH2 ; CH2  HH2 ; CH2  CZ3 ; CZ3  HZ3
CE3  CZ3 ; C    O   ; C    +R

-R   CA   N    H   ; CA   +R   C    O   ; CD1  CE2  NE1  HE1                                             
CH2  CE2  CZ2  HZ2 ; CZ2  CZ3  CH2  HH2 ; CE3  CH2  CZ3  HZ3                                             
CZ3  CD2  CE3  HE3 ; CG   HD1  CD1  NE1 ; CD1  CD2  CG   CB

Residue TYR
 21 23  8
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.005
HB1   HC     0.06
HB2   HC     0.06
CG    CA    -0.115
CD1   CA    -0.115
HD1   HA     0.115
CE1   CA    -0.115
HE1   HA     0.115
CZ    CA     0.150
OH    OHT   -0.585
HH    HO     0.435
CE2   CA    -0.115
HE2   HA     0.115
CD2   CA    -0.115
HD2   HA     0.115
C     C      0.50
O     O     -0.50
                                                                
-R   N   ; N    H   ; N    CA ; CA   HA  ; CA   CB  ; CA   C
CB   HB1 ; CB   HB2 ; CB   CG ; CG   CD1 ; CD1  HD1 ; CD1  CE1
CE1  HE1 ; CE1  CZ  ; CZ   OH ; OH   HH  ; CZ   CE2 ; CE2  HE2
CE2  CD2 ; CD2  HD2 ; CD2  CG ; C    O   ; C    +R

-R   CA   N    H   ; CA   +R   C    O   ; CE2  CG   CD2  HD2 ; CD2  CZ   CE2  HE2
CD1  CZ   CE1  HE1 ; CE1  CG   CD1  HD1 ; CD1  CD2  CG   CB  ; CE1  CE2  CZ   OH

Residue VAL
 16 17  2
N     N     -0.50
H     H      0.30
CA    CT     0.14
HA    HC     0.06
CB    CT    -0.06
HB    HC     0.06
CG1   CT    -0.18
HG11  HC     0.06
HG12  HC     0.06
HG13  HC     0.06
CG2   CT    -0.18
HG21  HC     0.06
HG22  HC     0.06
HG23  HC     0.06
C     C      0.50
O     O     -0.50

-R   N    ; N    H    ; N    CA   ; CA   HA   ; CA   CB   ; CA   C
CB   HB   ; CB   CG1  ; CB   CG2  ; CG1  HG11 ; CG1  HG12 ; CG1  HG13
CG2  HG21 ; CG2  HG22 ; CG2  HG23 ; C    O    ; C    +R

-R   CA   N    H ; CA   +R   C    O                                               

Residue CL 
  1  0  0
CL    CL    -1.0

Residue NA
  1  0  0
NA    SOD    1.0

Residue ZN
  1  0  0
ZN       ZN      +2.0

Residue Water
  3  3  0
OH2      OW      -0.834
H1       HW       0.417
H2       HW       0.417

OH2 H1 ; OH2 H2 ; H1 H2

!-------------------------------------------------------------------------------
"""

opls_bond = """
Variants

Variant XC_TERMINAL
  2  2  0  2  0

C ; O

C     C    0.5
O     O   -0.5

CA C ; C  O

Variant XN_TERMINAL
  2  2  0  2  0

N ; H

N     N   -0.50
H     H    0.30

N  H ; N  CA

Variant XN_TERMINAL PRO
  1  1  0  1  0

N

N     N   -0.14

N  CA

Variant C_TERMINAL
  2  3  1  3  1

C ; O

C     CCA   0.7
O     O2   -0.8
OXT   O2   -0.8

CA    0.04

CA C ; C  O ; C  OXT

CA  OXT  C    O                                               

Variant C_TERMINAL GLY
  2  3  1  3  1

C ; O

C     C     0.7
O     O2   -0.8
OXT   O2   -0.8

CA   -0.02

CA C ; C  O ; C  OXT

CA  OXT  C    O                                               

Variant N_TERMINAL
  2  4  1  4  0

N ; H

N     N3   -0.30
H1    H3    0.33
H2    H3    0.33
H3    H3    0.33

CA    0.25

N  H1 ; N  H2 ; N  H3 ; N  CA

Variant N_TERMINAL GLY
  2  4  1  4  0

N ; H

N     N3   -0.30
H1    H3    0.33
H2    H3    0.33
H3    H3    0.33

CA    0.19

N  H1 ; N  H2 ; N  H3 ; N  CA

Variant N_TERMINAL PRO
  1  3  5  3  0

N 

N     N    -0.07
H1    H     0.24
H2    H     0.24

CA    0.16
CD    0.16
HA    0.09 
HD1   0.09
HD2   0.09

N  H1 ; N  H2 ; N  CA

End

Links

Link DISULPHIDE_BRIDGE
  2  1  1  2  0

HG ; SG

SG    S    -0.2175

CB    0.0975

CB   SG ; SG   *R

  2  1  1  2  0

HG ; SG

SG    S    -0.2175

CB    0.0975

CB   SG ; SG   *R

End

Parameters

Bonds
C    OS     214.0    1.327
OS   CT     320.0    1.410
CT   CT     268.0    1.529
CT   HC     340.0    1.090
CT   CA     317.0    1.510
CT   OH     320.0    1.410
CT   SH     222.0    1.810
CT   S      222.0    1.810
CT   C      317.0    1.522
CT   N      337.0    1.449
CT   CCA    317.0    1.522
CT   N3     367.0    1.471
CT   N2     337.0    1.463
CT   C*     317.0    1.495
CT   CV     317.0    1.504
CT   CW     317.0    1.504
CT   CTP    268.0    1.529
HC   CTP    340.0    1.090
CA   CA     469.0    1.400
CA   HA     367.0    1.080
CA   OHT    450.0    1.364
CA   C*     546.0    1.352
CA   CB     469.0    1.400
CA   CN     469.0    1.400
CA   NA     427.0    1.381
HA   CP     367.0    1.080
HA   CV     367.0    1.080
HA   CW     367.0    1.080
OH   HO     553.0    0.945
HO   OHT    553.0    0.945
SH   HS     274.0    1.336
S    S      166.0    2.038
C    O      570.0    1.229
C    N      490.0    1.335
N    H      434.0    1.010
N    CTP    337.0    1.449
H    NA     434.0    1.010
CCA  O2     656.0    1.250
N3   H3     434.0    1.010
H3   N2     434.0    1.010
N2   CAA    481.0    1.340
C*   CB     388.0    1.459
CB   CN     447.0    1.419
CN   NA     428.0    1.380
NA   CP     477.0    1.343
NA   CW     427.0    1.381
CP   NB     488.0    1.335
CV   CW     518.0    1.371
CV   NB     410.0    1.394
C    OHC    450.0    1.364
HO   OHC    553.0    0.945
CW   CW     518.0    1.371
HW   OW     529.6    0.957
HW   HW      38.3    1.514
! ------------------------
"""

opls_angl = """
Angles
CT   C    OS     81.0    111.40
C    OS   CT     83.0    116.90
OS   CT   HC     35.0    109.50
O    C    OS     83.0    123.40
CT   CT   CT     58.4    112.70
CT   CT   HC     37.5    110.70
CT   CT   CA     63.0    114.00
CT   CT   OH     50.0    109.50
CT   CT   SH     50.0    108.60
CT   CT   S      50.0    114.70
CT   CT   C      63.0    111.10
CT   CT   N      80.0    109.70
CT   CT   CCA    63.0    111.10
CT   CT   N3     80.0    111.20
CT   CT   N2     80.0    111.20
CT   CT   C*     63.0    115.60
CT   CT   CV     63.0    114.00
CT   CT   CW     63.0    114.00
CT   CT   CTP    58.4    112.70
HC   CT   HC     33.0    107.80
HC   CT   CA     35.0    109.50
HC   CT   OH     35.0    109.50
HC   CT   SH     35.0    109.50
HC   CT   S      35.0    109.50
HC   CT   C      35.0    109.50
HC   CT   N      35.0    109.50
HC   CT   CCA    35.0    109.50
HC   CT   N3     35.0    109.50
HC   CT   N2     35.0    109.50
HC   CT   C*     35.0    109.50
HC   CT   CV     35.0    109.50
HC   CT   CW     35.0    109.50
HC   CT   CTP    37.5    110.70
C    CT   N      63.0    110.10
C    CT   N3     80.0    111.20
N    CT   CCA    63.0    110.10
CCA  CT   N3     80.0    111.20
CT   CA   CA     70.0    120.00
CA   CA   CA     63.0    120.00
CA   CA   HA     35.0    120.00
CA   CA   OHT    70.0    120.00
CA   CA   CB     63.0    120.00
CA   CA   CN     85.0    120.00
HA   CA   C*     35.0    120.00
HA   CA   CB     35.0    120.00
HA   CA   CN     35.0    120.00
HA   CA   NA     35.0    121.60
C*   CA   NA     70.0    108.70
CT   OH   HO     55.0    108.50
CA   OHT  HO     35.0    113.00
CT   SH   HS     44.0     96.00
CT   S    CT     62.0     98.90
CT   S    S      68.0    103.70
CT   C    O      80.0    120.40
CT   C    N      70.0    116.60
O    C    N      80.0    122.90
CT   N    CT     50.0    118.00
CT   N    C      50.0    121.90
CT   N    H      38.0    118.40
CT   N    CTP    50.0    118.00
C    N    H      35.0    119.80
C    N    CTP    50.0    121.90
H    N    H      35.0    120.00
CT   CCA  O2     70.0    117.00
O2   CCA  O2     80.0    126.00
CT   N3   CT     50.0    113.00
CT   N3   H3     35.0    109.50
H3   N3   H3     35.0    109.50
CT   N2   H3     35.0    118.40
CT   N2   CAA    50.0    123.20
H3   N2   H3     35.0    120.00
H3   N2   CAA    35.0    120.00
N2   CAA  N2     70.0    120.00
CT   C*   CA     70.0    125.00
CT   C*   CB     70.0    128.60
CA   C*   CB     85.0    106.40
CA   CB   C*     85.0    134.90
CA   CB   CN     85.0    116.20
C*   CB   CN     85.0    108.80
CA   CN   CB     85.0    122.70
CA   CN   NA     70.0    132.80
CB   CN   NA     70.0    104.40
CA   NA   H      35.0    120.00
CA   NA   CN     70.0    111.60
H    NA   CN     35.0    123.10
H    NA   CP     35.0    126.35
H    NA   CW     35.0    126.35
CP   NA   CW     70.0    107.30
HA   CP   NA     35.0    120.00
HA   CP   NB     35.0    120.00
NA   CP   NA     70.0    110.75
NA   CP   NB     70.0    111.60
CT   CV   CW     70.0    130.70
HA   CV   CW     35.0    128.20
CT   CV   NB     70.0    124.50
HA   CV   NB     35.0    120.00
CW   CV   NB     70.0    108.30
CT   CW   NA     70.0    121.60
CT   CW   CV     70.0    130.70
CT   CW   CW     70.0    130.70
HA   CW   NA     35.0    121.60
HA   CW   CV     35.0    130.70
HA   CW   CW     35.0    130.70
NA   CW   CV     70.0    108.70
NA   CW   CW     70.0    106.30
CP   NB   CV     70.0    105.30
CT   CTP  HC     37.5    110.70
CT   CTP  N      80.0    109.70
HC   CTP  HC     33.0    107.80
HC   CTP  N      35.0    109.50
CT   C    OHC    80.0    120.40  
O    C    OHC    80.0    121.00
C    OHC  HO     35.0    113.50
HW   OW   HW     34.1    104.52
HW   HW   OW      0.0     37.74
! -----------------------------
"""

opls_dihe = """
Dihedrals
OS   C    CT   N       0.000    0.000    0.820    0.000
CT   C    OS   CT      0.000    0.000    4.830    0.000
O    C    OS   CT      0.000    0.000    4.830    0.000
OS   C    CT   HC      0.000    0.000    0.000    0.000
HC   CT   OS   C       0.000    0.000    0.000    0.000
OS   C    CT   CT      0.000    0.000    0.820    0.000
CT   CT   CT   CT      0.000    1.740   -0.157    0.279
CT   CT   CT   HC      0.000    0.000    0.000    0.366
CT   CT   CT   S       0.000    2.619   -0.620    0.258
CT   CT   CT   C       0.000   -1.697   -0.456    0.585
CT   CT   CT   N       0.000    0.845   -0.962    0.713
CT   CT   CT   CCA     0.000   -3.185   -0.825    0.493
CW   CT   CT   CCA     0.000   -3.185   -0.825    0.493
CT   CT   CT   N3      0.000    2.732   -0.229    0.485
CT   CT   CT   N2      0.000    1.964    0.000    0.659
CT   CT   CT   CTP     0.000    1.740   -0.157    0.279
HC   CT   CT   HC      0.000    0.000    0.000    0.318
HC   CT   CT   CA      0.000    0.000    0.000    0.462
HC   CT   CT   OH      0.000    0.000    0.000    0.468
HC   CT   CT   SH      0.000    0.000    0.000    0.452
HC   CT   CT   S       0.000    0.000    0.000    0.452
HC   CT   CT   C       0.000    0.000    0.000   -0.076
HC   CT   CT   N       0.000    0.000    0.000    0.464
HC   CT   CT   CCA     0.000    0.000    0.000   -0.225
HC   CT   CT   N3      0.000    0.000    0.000    0.384
HC   CT   CT   N2      0.000    0.000    0.000    0.464
HC   CT   CT   C*      0.000    0.000    0.000    0.462
HC   CT   CT   CV      0.000    0.000    0.000    0.462
HC   CT   CT   CW      0.000    0.000    0.000    0.462
HC   CT   CT   CTP     0.000    0.000    0.000    0.366
CA   CT   CT   C       0.000   -1.697   -0.456    0.585
CA   CT   CT   N       0.000    0.845   -0.962    0.713
CA   CT   CT   CCA     0.000   -1.697   -0.456    0.585
CA   CT   CT   N3      0.000    0.845   -0.962    0.713
OH   CT   CT   C       0.000   -6.180    0.000    0.000
OH   CT   CT   N       0.000    6.280   -1.467    2.030
OH   CT   CT   CCA     0.000   -6.180    0.000    0.000
OH   CT   CT   N3      0.000    6.280   -1.467    2.030
SH   CT   CT   C       0.000   -4.214   -2.114    0.969
SH   CT   CT   N       0.000    0.583   -1.163    0.141
SH   CT   CT   CCA     0.000   -4.214   -2.114    0.969
SH   CT   CT   N3      0.000    0.583   -1.163    0.141
S    CT   CT   C       0.000   -4.214   -2.114    0.969
S    CT   CT   N       0.000    0.583   -1.163    0.141
S    CT   CT   CCA     0.000   -4.214   -2.114    0.969
S    CT   CT   N3      0.000    0.583   -1.163    0.141
C    CT   CT   C       0.000   -1.697   -0.456    0.585
C    CT   CT   N       0.000    0.845   -0.962    0.713
C    CT   CT   CCA     0.000   -1.697   -0.456    0.585
C    CT   CT   N3      0.000    0.845   -0.962    0.713
C    CT   CT   C*      0.000   -1.697   -0.456    0.585
C    CT   CT   CV      0.000   -1.697   -0.456    0.585
C    CT   CT   CW      0.000   -1.697   -0.456    0.585
N    CT   CT   CCA     0.000    0.845   -0.962    0.713
N    CT   CT   C*      0.000    0.845   -0.962    0.713
N    CT   CT   CV      0.000    0.845   -0.962    0.713
N    CT   CT   CW      0.000    0.845   -0.962    0.713
CCA  CT   CT   CCA     0.000    0.845   -0.962    0.713
CCA  CT   CT   N3      0.000    0.845   -0.962    0.713
CCA  CT   CT   C*      0.000   -1.697   -0.456    0.585
N3   CT   CT   C*      0.000    0.845   -0.962    0.713
CT   CT   CA   CA      0.000    0.000    0.000    0.000
HC   CT   CA   CA      0.000    0.000    0.000    0.000
CT   CT   OH   HO      0.000   -0.356   -0.174    0.492
HC   CT   OH   HO      0.000    0.000    0.000    0.450
CT   CT   SH   HS      0.000   -0.759   -0.282    0.603
HC   CT   SH   HS      0.000    0.000    0.000    0.451
CT   CT   S    CT      0.000    0.925   -0.576    0.677
CT   CT   S    S       0.000    1.941   -0.836    0.935
HC   CT   S    CT      0.000    0.000    0.000    0.647
HC   CT   S    S       0.000    0.000    0.000    0.558
CT   CT   C    O       0.000    0.000    0.000    0.000
CT   CT   C    N       0.000    1.173    0.189   -1.200
HC   CT   C    O       0.000    0.000    0.000    0.000
HC   CT   C    N       0.000    0.000    0.000    0.000
N    CT   C    O       0.000    0.000    0.000    0.000
N    CT   C    N       0.000    1.816    1.222    1.581
N3   CT   C    O       0.000    0.000    0.000    0.000
N3   CT   C    N       0.000    1.816    1.222    1.581
CT   CT   N    C       0.000    0.000    0.462    0.000
CT   CT   N    H       0.000    0.000    0.000    0.000
CT   CT   N    CTP     0.000    2.377   -0.367    0.000
HC   CT   N    CT      0.000    0.000    0.000    0.000
HC   CT   N    C       0.000    0.000    0.000    0.000
HC   CT   N    H       0.000    0.000    0.000    0.000
HC   CT   N    CTP     0.000    0.000    0.000    0.000
C    CT   N    C       0.000   -2.365    0.912   -0.850
C    CT   N    H       0.000    0.000    0.000    0.000
C    CT   N    CTP     0.000   -0.869    0.626   -1.751
CCA  CT   N    C       0.000   -2.365    0.912   -0.850
CCA  CT   N    H       0.000    0.000    0.000    0.000
CCA  CT   N    CTP     0.000   -0.869    0.626   -1.751
CT   CT   CCA  O2      0.000    0.000    0.820    0.000
HC   CT   CCA  O2      0.000    0.000    0.000    0.000
N    CT   CCA  O2      0.000    0.000    0.000    0.000
N3   CT   CCA  O2      0.000    0.000    0.000    0.000
CT   CT   N3   CT      0.000    1.740   -0.157    0.279
CT   CT   N3   H3      0.000    0.000    0.000    0.347
HC   CT   N3   CT      0.000    0.000    0.000    0.366
HC   CT   N3   H3      0.000    0.000    0.000    0.261
C    CT   N3   CT      0.000   -1.697   -0.456    0.585
C    CT   N3   H3      0.000    0.000    0.000    0.347
CCA  CT   N3   H3      0.000    0.000    0.000    0.347
CT   CT   N2   H3      0.000    0.000    0.000    0.000
CT   CT   N2   CAA     0.000    1.829    0.243   -0.498
HC   CT   N2   H3      0.000    0.000    0.000    0.000
HC   CT   N2   CAA     0.000    0.000    0.000    0.177
CT   CT   C*   CA      0.000   -0.714    0.000    0.000
CT   CT   C*   CB      0.000    0.000    0.000    0.000
HC   CT   C*   CA      0.000    0.000    0.000   -0.480
HC   CT   C*   CB      0.000    0.000    0.000    0.000
CT   CT   CV   CW      0.000    0.000    0.000    0.000
CT   CT   CV   NB      0.000    2.366   -0.262    0.505
HC   CT   CV   CW      0.000    0.000    0.000    0.000
HC   CT   CV   NB      0.000    0.000    0.000    0.419
CT   CT   CW   NA      0.000    2.366   -0.262    0.505
CT   CT   CW   CV      0.000    0.000    0.000    0.000
CT   CT   CW   CW      0.000    0.000    0.000    0.000
HC   CT   CW   NA      0.000    0.000    0.000    0.419
HC   CT   CW   CV      0.000    0.000    0.000    0.000
HC   CT   CW   CW      0.000    0.000    0.000    0.000
CT   CT   CTP  HC      0.000    0.000    0.000    0.366
CT   CT   CTP  N       0.000    0.845   -0.962    0.713
HC   CT   CTP  HC      0.000    0.000    0.000    0.318
HC   CT   CTP  N       0.000    0.000    0.000    0.464
CT   CA   CA   CA      0.000    0.000    3.625    0.000
CT   CA   CA   HA      0.000    0.000    3.625    0.000
CA   CA   CA   CA      0.000    0.000    3.625    0.000
CA   CA   CA   HA      0.000    0.000    3.625    0.000
CA   CA   CA   OHT     0.000    0.000    3.625    0.000
CA   CA   CA   CB      0.000    0.000    3.625    0.000
CA   CA   CA   CN      0.000    0.000    3.625    0.000
HA   CA   CA   HA      0.000    0.000    3.625    0.000
HA   CA   CA   OHT     0.000    0.000    3.625    0.000
HA   CA   CA   CB      0.000    0.000    3.625    0.000
HA   CA   CA   CN      0.000    0.000    3.625    0.000
CA   CA   OHT  HO      0.000    0.000    1.682    0.000
HA   CA   C*   CT      0.000    0.000    6.525    0.000
HA   CA   C*   CB      0.000    0.000    6.525    0.000
NA   CA   C*   CT      0.000    0.000    6.525    0.000
NA   CA   C*   CB      0.000    0.000    6.525    0.000
CA   CA   CB   C*      0.000    0.000    3.500    0.000
CA   CA   CB   CN      0.000    0.000    3.500    0.000
CA   CA   CN   CB      0.000    0.000    3.625    0.000
CA   CA   CN   NA      0.000    0.000    3.625    0.000
HA   CA   CB   C*      0.000    0.000    3.500    0.000
HA   CA   CB   CN      0.000    0.000    3.500    0.000
HA   CA   CN   CB      0.000    0.000    3.625    0.000
HA   CA   CN   NA      0.000    0.000    3.625    0.000
HA   CA   NA   H       0.000    0.000    1.500    0.000
HA   CA   NA   CN      0.000    0.000    1.500    0.000
C*   CA   NA   H       0.000    0.000    1.500    0.000
C*   CA   NA   CN      0.000    0.000    1.500    0.000
CT   S    S    CT      0.000    0.000   -7.414    1.705
CT   C    N    CT      0.000    2.800    6.089    0.000
CT   C    N    H       0.000    0.000    4.900    0.000
CT   C    N    CTP     0.000    2.800    6.089    0.000
O    C    N    CT      0.000    0.000    6.089    0.000
O    C    N    H       0.000    0.000    4.900    0.000
O    C    N    CTP     0.000    0.000    6.089    0.000
CT   N    CTP  CT      0.000    1.430    1.029   -5.633
CT   N    CTP  HC      0.000    0.000    0.000    0.000
C    N    CTP  CT      0.000    0.000    0.462    0.000
C    N    CTP  HC      0.000    0.000    0.000    0.000
CT   N2   CAA  N2      0.000    0.000    7.936    0.000
H3   N2   CAA  N2      0.000    0.000    3.900    0.000
CT   C*   CB   CA      0.000    0.000    1.675    0.000
CT   C*   CB   CN      0.000    0.000    1.675    0.000
CA   C*   CB   CA      0.000    0.000    1.675    0.000
CA   C*   CB   CN      0.000    0.000    1.675    0.000
CA   CB   CN   CA      0.000    0.000    3.000    0.000
CA   CB   CN   NA      0.000    0.000    3.000    0.000
C*   CB   CN   CA      0.000    0.000    3.000    0.000
C*   CB   CN   NA      0.000    0.000    3.000    0.000
CA   CN   NA   CA      0.000    0.000    1.525    0.000
CA   CN   NA   H       0.000    0.000    1.525    0.000
CB   CN   NA   CA      0.000    0.000    1.525    0.000
CB   CN   NA   H       0.000    0.000    1.525    0.000
H    NA   CP   HA      0.000    0.000    2.325    0.000
H    NA   CP   NA      0.000    0.000    2.325    0.000
H    NA   CP   NB      0.000    0.000    2.325    0.000
CW   NA   CP   HA      0.000    0.000    2.325    0.000
CW   NA   CP   NA      0.000    0.000    2.325    0.000
CW   NA   CP   NB      0.000    0.000    2.325    0.000
H    NA   CW   CT      0.000    0.000    1.400    0.000
H    NA   CW   HA      0.000    0.000    1.600    0.000
H    NA   CW   CV      0.000    0.000    1.400    0.000
H    NA   CW   CW      0.000    0.000    1.400    0.000
CP   NA   CW   CT      0.000    0.000    1.400    0.000
CP   NA   CW   HA      0.000    0.000    1.600    0.000
CP   NA   CW   CV      0.000    0.000    1.400    0.000
CP   NA   CW   CW      0.000    0.000    1.400    0.000
HA   CP   NB   CV      0.000    0.000    5.000    0.000
NA   CP   NB   CV      0.000    0.000    5.000    0.000
CT   CV   CW   HA      0.000    0.000    5.375    0.000
CT   CV   CW   NA      0.000    0.000    5.375    0.000
HA   CV   CW   CT      0.000    0.000    5.375    0.000
HA   CV   CW   NA      0.000    0.000    5.375    0.000
NB   CV   CW   CT      0.000    0.000    5.375    0.000
NB   CV   CW   HA      0.000    0.000    5.375    0.000
NB   CV   CW   NA      0.000    0.000    5.375    0.000
CT   CV   NB   CP      0.000    0.000    2.400    0.000
HA   CV   NB   CP      0.000    0.000    2.400    0.000
CW   CV   NB   CP      0.000    0.000    2.400    0.000
CT   CW   CW   HA      0.000    0.000    5.375    0.000
CT   CW   CW   NA      0.000    0.000    5.375    0.000
HA   CW   CW   NA      0.000    0.000    5.375    0.000
NA   CW   CW   NA      0.000    0.000    5.375    0.000
OHC  C    CT   CT      0.000    0.000    0.820    0.000
OHC  C    CT   HC      0.000    0.000    0.000    0.000
O    C    OHC  HO      0.000    0.000    4.830    0.000
CT   C    OHC  HO      0.000    0.000    4.830    0.000
! -----------------------------------------------------
"""

opls_impr = """
Impropers
O    C    OS   CT      0.000    0.000   10.500    0.000
CA   CA   CA   CT      0.000    0.000    1.100    0.000
CA   CA   CA   HA      0.000    0.000    1.100    0.000
CA   CA   CA   OHT     0.000    0.000    1.100    0.000
CA   CN   CA   HA      0.000    0.000    1.100    0.000
CA   CB   CA   HA      0.000    0.000    1.100    0.000
C*   HA   CA   NA      0.000    0.000    1.100    0.000
CT   N    C    O       0.000    0.000   10.500    0.000
C    CT   N    H       0.000    0.000    1.000    0.000
C    CT   N    CTP     0.000    0.000    1.000    0.000
C    H    N    H       0.000    0.000    1.000    0.000
CT   O2   CCA  O2      0.000    0.000   10.500    0.000
CA   CB   C*   CT      0.000    0.000    1.100    0.000
C*   CA   CB   CN      0.000    0.000    1.100    0.000
CA   CB   CN   NA      0.000    0.000    1.100    0.000
CA   CN   NA   H       0.000    0.000    1.000    0.000
CP   CW   NA   H       0.000    0.000    1.000    0.000
HA   NA   CP   NA      0.000    0.000    1.100    0.000
HA   NA   CP   NB      0.000    0.000    1.100    0.000
CT   CW   CV   NB      0.000    0.000    1.100    0.000
CW   HA   CV   NB      0.000    0.000    1.100    0.000
CT   CV   CW   NA      0.000    0.000    1.100    0.000
CV   HA   CW   NA      0.000    0.000    1.100    0.000
CT   CW   CW   NA      0.000    0.000    1.100    0.000
CW   HA   CW   NA      0.000    0.000    1.100    0.000
N2   N2   CAA  N2      0.000    0.000   10.500    0.000
CAA  CT   N2   H3      0.000    0.000    1.000    0.000
CAA  H3   N2   H3      0.000    0.000    1.000    0.000
CT   O    C    OHC     0.000    0.000   10.500    0.000
! -----------------------------------------------------
"""


mass = { 1 : 1.00794, 2 : 4.00260, 3 : 6.94100, 4 : 9.01218, 5 : 10.8110, 6 : 12.0107, 7 : 14.0067, 8 : 15.9994, 9 : 18.9984, 10 : 20.1797,
	11 : 22.9898, 12 : 24.3050, 13 : 26.9815, 14 : 28.0855, 15 : 30.9738, 16 : 32.0650, 17 : 35.4530, 18 : 39.9480, 19 : 39.0983, 20 : 40.0780,
	21 : 44.9559, 22 : 47.8670, 23 : 50.9415, 24 : 51.9961, 25 : 54.9380, 26 : 55.8450, 27 : 58.9332, 28 : 58.6934, 29 : 63.5460, 30 : 65.3900,
	31 : 69.7230, 32 : 72.6400, 33 : 74.9216, 34 : 78.9600, 35 : 79.9040, 36 : 83.8000, 37 : 85.4678, 38 : 87.6200, 39 : 88.9059, 40 : 91.2240,
	41 : 92.9064, 42 : 95.9400, 43 : 98.9063, 44 : 101.0700, 45 : 102.9060, 46 : 106.4200, 47 : 107.8680, 48 : 112.4110, 49 : 114.8180, 50 : 118.7100,
	51 : 121.7600, 52 : 127.6000, 53 : 126.9040, 54 : 131.2930, 55 : 132.9050, 56 : 137.2370, 57 : 138.9050, 58 : 140.1160, 59 : 140.9080, 60 : 144.2400,
	61 : 146.9150, 62 : 150.3600, 63 : 151.9640, 64 : 157.2500, 65 : 158.9250, 66 : 162.5000, 67 : 164.9300, 68 : 167.2590, 69 : 168.9340, 70 : 173.0400,
	71 : 174.9670, 72 : 178.4900, 73 : 180.9480, 74 : 183.8400, 75 : 186.2070, 76 : 190.2300, 77 : 192.2170, 78 : 195.0780, 79 : 196.9670, 80 : 200.5900,
	81 : 204.3830, 82 : 207.2000, 83 : 208.9800, 84 : 208.9820, 85 : 209.9870, 86 : 222.0180, 87 : 223.0200, 88 : 226.0250, 89 : 227.0280, 90 : 232.0380,
	91 : 231.0360, 92 : 238.0290, 93 : 237.0480, 94 : 244.0640, 95 : 243.0610, 96 : 247.0700, 97 : 247.0700, 98 : 251.0800, 99 : 252.0830, 100 : 257.0950,
	101 : 258.0990, 102 : 259.1010, 103 : 262.1100, 104 : 261.1090, 105 : 262.1140, 106 : 263.1190, 107 : 262.1230, 108 : 265.1310, 109 : 266.1380 }


atom_type = collections.OrderedDict()
resi_name = None
resi_atom = []
resi_type = []
resi_chrg = []
resi_bond = []
resi_impr = []
f = open( sys.argv[1], "rt" )
l = f.readline()
while( l != "" ):
	t = l.upper().strip().split()
	if( len( t ) == 4 and t[0] == "MASS" ):
		m = float( t[3] )
		atom_type[t[2]] = sorted( [ ( math.fabs( m - mass[i] ), i ) for i in list( mass ) ] )[0][1]
	elif( len( t ) == 3 and t[0] == "RESI" ):
		resi_name = t[1]
	elif( len( t ) == 4 and t[0] == "ATOM" ):
		resi_atom.append( t[1] )
		resi_type.append( t[2] )
		resi_chrg.append( float( t[3] ) )
	elif( len( t ) >= 3 and t[0] == "BOND" ):
		resi_bond.append( " ".join( t[1:3] ) )
	elif( len( t ) >= 5 and t[0] == "IMPH" ):
		resi_impr.append( " ".join( t[1:5] ) ) #( t[1], t[2], t[3], t[4] ) )
	l = f.readline()
f.close()


ffld = { "bond": collections.OrderedDict(),
		"angl": collections.OrderedDict(),
		"dihe": collections.OrderedDict(),
		"impr": collections.OrderedDict(),
		"nbnd": collections.OrderedDict() }
f = open( sys.argv[2], "rt" )
key = None
pat = re.compile( "([^\ ]+)[\ ]+[0-9\.]+[\ ]+([0-9\.\-]+)[\ ]+([0-9\.]+)" )
cte = 2. / math.pow( 2., 1. / 6. )
l = f.readline()
while( l != "" ):
	if( l[0:4].upper() == "BOND" ):
		key = "bond"
	elif( l[0:4].upper() == "ANGL" ):
		key = "angl"
	elif( l[0:4].upper() == "DIHE" ):
		key = "dihe"
	elif( l[0:4].upper() in [ "IMPH", "IMPR" ] ):
		key = "impr"
	elif( l[0:4].upper() == "NONB" ):
		key = "nbnd"
	# --------------------------------------
	t = l.upper().strip().split()
	if( key == "bond" and len( t ) >= 4 ):
		ffld[key][" ".join( t[0:2] )] = ( float( t[2] ), float( t[3] ) )
	elif( key == "angl" and len( t ) >= 5 ):
		ffld[key][" ".join( t[0:3] )] = ( float( t[3] ), float( t[4] ) )
	elif( key in [ "dihe", "impr" ] and len( t ) >= 7 ):
		k = " ".join( t[0:4] )
		if( ffld[key].has_key( k ) ):
			ffld[key][k] += [ float( i ) for i in t[4:7] ]
		else:
			ffld[key][k] = [ float( i ) for i in t[4:7] ]
	elif( key == "nbnd" and len( t ) >= 3 ):
		t = pat.findall( l.upper() )
		if( len( t ) > 0 ):
			ffld[key][t[0][0]] = ( cte * float( t[0][2] ), math.fabs( float( t[0][1] ) ) )
	l = f.readline()
f.close()



def __dihe( v, x, y ):
	f = 0
	g = [ .0, .0, .0, .0 ]
	n = len( x )
	for i in range( n ):
		t_1   = ( 1.0 + math.cos( x[i] ) )
		t_2   = ( 1.0 - math.cos( 2.0 * x[i] ) )
		t_3   = ( 1.0 + math.cos( 3.0 * x[i] ) )
		t     = 0.5 * ( v[0] + v[1] * t_1 + v[2] * t_2 + v[3] * t_3 )
		d     = ( t - y[i] )
		f    += d * d
		g[0] += d
		g[1] += d * t_1
		g[2] += d * t_2
		g[3] += d * t_3
	return( f, g )


def __impr( v, x, y ):
	f = 0
	g = [ .0, .0, .0, .0 ]
	n = len( x )
	for i in range( n ):
		t_1   = ( 1.0 + math.cos( x[i] ) )
		t_2   = ( 1.0 - math.cos( 2.0 * x[i] ) )
		t     = 0.5 * ( v[0] + v[1] * t_1 + v[2] * t_2 )
		d     = ( t - y[i] )
		f    += d * d
		g[0] += d
		g[1] += d * t_1
		g[2] += d * t_2
	return( f, g )


def __fire( function, x, y ):
	msiz = 0.1
	coor = [ .0, .0, .0, .0 ]
	nstp = 0
	alph = 0.1
	ssiz = msiz
	velo = [ .0, .0, .0, .0 ]
	step = [ .0, .0, .0, .0 ]
	func, grad = function( coor, x, y )
	norm = math.sqrt( grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2] + grad[3] * grad[3] )
	it   = 0
	while( it < 5000 and norm > 1.e-8 ):
		vsiz = math.sqrt( velo[0] * velo[0] + velo[1] * velo[1] + velo[2] * velo[2] + velo[3] * velo[3] )
		vfac = - ( velo[0] * grad[0] + velo[1] * grad[1] + velo[2] * grad[2] + velo[3] * grad[3] )
		if( vfac > 0.0 ):
			for i in [0, 1, 2, 3]:
				velo[i] = ( 1.0 - alph ) * velo[i] - alph * grad[i] / norm * vsiz
			if( nstp > 5 ):
				ssiz  = min( ssiz * 1.1, msiz )
				alph *= 0.99
			nstp += 1
		else:
			velo  = [ .0, .0, .0, .0 ]
			alph  = 0.1
			nstp  = 0
			ssiz *= 0.5
		for i in [0, 1, 2, 3]:
			velo[i] -= ssiz * grad[i]
			step[i]  = ssiz * velo[i]
		disp = math.sqrt( step[0] * step[0] + step[1] * step[1] + step[2] * step[2] + step[3] * step[3] )
		if( disp > ssiz ):
			for i in [0, 1, 2, 3]:
				step[i] *= ssiz / disp
		for i in [0, 1, 2, 3]:
			coor[i] += step[i]
		func, grad = function( coor, x, y )
		norm = math.sqrt( grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2] + grad[3] * grad[3] )
		it  += 1
	print( "%20.10lf%20.10lf"%( func, norm ) )
	return( coor, func, norm )
	



N  = 100
_d = 2. * math.pi / float( N - 1 )
X  = [ - math.pi + float( i ) * _d for i in range( N ) ]

f = open( "opls", "wt" )
f.write( opls_type )
for k in list( atom_type ):
	f.write( "x%-4s%4d%10.5lf%10.5lf\n"%( k, atom_type[k], ffld["nbnd"][k][0], ffld["nbnd"][k][1] ) )
f.write( "End\n" )
f.write( opls_resi )
f.write( "Residue %s\n"%( resi_name ) )
n = len( resi_atom )
f.write( "%d %d %d\n"%( n, len( resi_bond ), len( resi_impr ) ) )
for i in range( n ):
	f.write( "%-9sx%-5s%8.3lf\n"%( resi_atom[i], resi_type[i], resi_chrg[i] ) )
f.write( "\n" )
f.write( "\n".join( resi_bond ) )
f.write( "\n\n" )
f.write( "\n".join( resi_impr ) )
f.write( "\nEnd\n" )
f.write( opls_bond )
for k in list( ffld["bond"] ):
	t = k.split()
	f.write( "x%-4sx%-4s%7.1lf%9.3lf\n"%( t[0], t[1], ffld["bond"][k][0], ffld["bond"][k][1] ) )
f.write( "End\n" )
f.write( opls_angl )
for k in list( ffld["angl"] ):
	t = k.split()
	f.write( "x%-4sx%-4sx%-4s%6.1lf%10.2lf\n"%( t[0], t[1], t[2], ffld["angl"][k][0], ffld["angl"][k][1] ) )
f.write( "End\n" )
f.write( opls_dihe )
for k in list( ffld["dihe"] ):
	y = [ 0.0 for i in range( N ) ]
	for i in range( 0, len( ffld["dihe"][k] ), 3 ):
		d = ffld["dihe"][k][i+2] / 180.0 * math.pi
		for j in range( N ):
			y[j] += ffld["dihe"][k][i] * ( 1 + math.cos( ffld["dihe"][k][i+1] * X[j] + d ) )
	rx = []
	ry = []
	for i in range( N ):
		if( math.fabs( y[i] ) <= 200.0 ):
			rx.append( X[i] )
			ry.append( y[i] )
	print( "%-30s"%( k ), end = "" )
	coor, func, norm = __fire( __dihe, rx, ry )
	if( func >= 10.0 ):
		f.write( "!-[fit] func/norm: %20.10lf / %20.10lf\n"%( func, norm ) )
	t = k.split()
	for i in [0, 1, 2, 3]:
		if( t[i] != "X" ):
			f.write( "x%-4s"%( t[i] ) )
		else:
			f.write( "X    " )
	f.write( "%8.3lf%9.3lf%9.3lf%9.3lf\n"%( coor[0], coor[1], coor[2], coor[3] ) )	
f.write( "End\n" )
f.write( opls_impr )
for k in list( ffld["impr"] ):
	y = [ 0.0 for i in range( N ) ]
	for i in range( 0, len( ffld["impr"][k] ), 3 ):
		d = ffld["impr"][k][i+2] / 180.0 * math.pi
		for j in range( N ):
			y[j] += ffld["impr"][k][i] * ( 1 + math.cos( ffld["impr"][k][i+1] * X[j] + d ) )
	rx = []
	ry = []
	for i in range( N ):
		if( math.fabs( y[i] ) <= 200.0 ):
			rx.append( X[i] )
			ry.append( y[i] )
	print( "%-30s"%( k ), end = "" )
	coor, func, norm = __fire( __impr, rx, ry )
	if( func >= 10.0 ):
		f.write( "!-[fit] func/norm: %20.10lf / %20.10lf\n"%( func, norm ) )
	t = k.split()
	for i in [0, 1, 2, 3]:
		if( t[i] != "X" ):
			f.write( "x%-4s"%( t[i] ) )
		else:
			f.write( "X    " )
	f.write( "%8.3lf%9.3lf%9.3lf%9.3lf\n"%( coor[0], coor[1], coor[2], coor[3] ) )	
f.write( "End\n\nEnd\nEnd" )
f.close()
