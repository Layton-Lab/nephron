# The values here are updated from the male human model. 

N=200
# physical constants
RT    = 2.57
RTosm = 1.93e4
F     = 96.5
visc  = 6.4e-6
# number of variables for H-K-ATPase
Natp = 14

# Vwbar = Molar volume of water  (cm3/mmole)
# Formula weight of H20 = 2*1.00797 + 15.9994 = 18.0153
# There is 1 mole H20/18.0153 gm of H20.
# Density of water at 37 C = 0.99335 gm/cm3 (CRC Handbook)
# Vwbar = (gm/millimole)/(gm/cm3)=(cm3/millimole)

Vwbar = 18.0153/1000/0.99335

# normalization constants
Cref = 1e-3  # 1 mM = 0.01 mmol/cm3
Vref = 1e-4  # in cm3/cm2 epith
Pfref = 1.0   # in cm/s
href = 1e-5  # in cm/s
EPref= 1e-3  # in volts

# acid-base pair balance constants
pKHCO3 = 3.57
pKHPO4 = 6.80
pKNH3 = 9.15
pKbuf = 7.5
pKHCO2 = 3.76

# transporter parameters for NKCC2 F isoform
poppnkccF = 3.9280e4
popnkccF = 3.578e5
pnkccpF = 1e4
pnkccppF = 1.098e3
pnmccpF = 2.0e3
pnmccppF = 219.6
bn2F = 58.93
bc2F = 13.12
bk2F = 9.149
bm2F = 9.149

# transporter parameters for NKCC2 A isoform
poppnkccA = 7.535e4
popnkccA = 2.594e5
pnkccpA = 1e4
pnkccppA = 2.904e3
pnmccpA = 2.0e3
pnmccppA = 580.8
bn2A = 118.8
bc2A = 0.08834
bk2A = 1.8710e4
bm2A = 1.8710e4

# transporter parameters for NKCC2 B isoform
poppnkccB = 2.517e5
popnkccB = 2.596e5
pnkccpB = 1e4
pnkccppB = 9.695e3
pnmccpB = 2.0e3
pnmccppB = 1.939e3
bn2B = 275.0
bc2B = 0.08157
bk2B = 5.577e3
bm2B = 5.577e3

# transporter parameters for KCC (KCC4 isoform)
poppkcc = 3.928e4
popkcc = 3.577e5
pkccp = 1.0e4
pkccpp = 1.098e3
pmccp = 2.0e3
pmccpp = 219.6
bckcc = 21.08
bkkcc = 1.45
bmkcc = 1.45

# transporter parameters for NCC
popncc = 4.295e6
poppncc = 0.1e6
pnpncc = 7692.0
pnppncc = 179.0
dKnncc = 0.293
dKcncc = 112.7
dKncncc = 0.565
dKcnncc = 1.47e-3

# transporter parameters for Pendrin
Pclepd = 10000.0
Pcliclepd = 1.239
Pbieclepd = 10.76
Poheclepd = 0.262
dKclpd = 3.01
dKbipd = 5.94
dKohpd = 1.38e-6

# transporter parameters for AE1
Pbp = 1247.0
Pbpp = 135.0
Pcp = 562.0
Pcpp = 61.0
dKbp = 198.0
dKbpp = 198.0
dKcp = 50.0
dKcpp = 50.0

# parameters
pHplasma = 7.4
phpap = 7.3
pKHCO3 = 3.57
pKHPO4 = 6.8
pKNH3 = 9.15
pKHCO2 = 3.76

# parameters used for NHE3i, NKCC2i, SNB, pregnancy

TotSodCM = 144.0
TotSodOI_noinhib = 299.0
TotSodOI_100NKCCinhib = 180.0
TotSodOI_70NKCCinhib = 239.5
TotSodOI_MP = 299.0*0.95
TotSodOI_LP = 299.0*0.95
TotSodPap_noinhib = 349.0
TotSodPap_50inhib = 324.0
TotSodPap_80inhib = 274.0
TotSodPap_100NKCCinhib = 140.0
TotSodPap_70NKCCinhib = 244.5
TotSodPap_MP = 349.0*0.95
TotSodPap_LP = 349.0*0.95

TotPotCM = 4.9
TotPotOI_noinhib = 10.0
TotPotOI_100NKCCinhib = 5.0
TotPotOI_70NKCCinhib = 7.5
TotPotOI_MP = 10.0*1.0
TotPotOI_LP = 10.0*1.15
TotPotPap_noinhib = 20.0
TotPotPap_80inhib = 10.0
TotPotPap_100NKCCinhib = 5.0
TotPotPap_70NKCCinhib = 12.5
TotPotPap_MP = 20.0*1.0
TotPotPap_LP = 20.0*1.15

TotCloCM = 116.17
TotCloOI_noinhib = 279.94
TotCloOI_100NKCCinhib = 158.0
TotCloOI_70NKCCinhib = 218.97
TotCloOI_MP = 279.94*0.95
TotCloOI_LP = 279.94*0.95
TotCloPap_noinhib = 344.92
TotCloPap_50inhib = 320.0
TotCloPap_80inhib = 260.0
TotCloPap_100NKCCinhib = 136.0
TotCloPap_70NKCCinhib = 240.46
TotCloPap_MP = 344.92*0.95
TotCloPap_LP = 344.92*0.95

TotureaCM = 8.0
TotureaOI_noinhib = 60.0
TotureaOI_50inhib = 50.0
TotureaOI_80inhib = 20.0
TotureaOI_100NKCCinhib = 20.0
TotureaOI_70NKCCinhib = 40.0

TotureaPap_noinhib = 200.0
TotureaPap_50inhib = 180.0
TotureaPap_80inhib = 20.0
TotureaPap_100NKCCinhib = 20.0
TotureaPap_70NKCCinhib = 110.0


TotHCO3OI_noinhib = 25.0
TotHCO3OI_MP = 25.0*0.95
TotHCO3OI_LP = 25.0*0.95
TotHCO3Pap_noinhib = 25.0
TotHCO3Pap_MP = 25.0*0.95
TotHCO3Pap_LP = 25.0*0.9645