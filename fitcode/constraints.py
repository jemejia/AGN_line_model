from fitcode.line import kms_to_wl 

#----------------- continuum limits -----------------
FUV_limits=[1216.0,2000.0]
UV_limits=[2000.0,4500.0]
OP_limits=[4000.0,8800.0]
#----------------- continuum limits ----------------





em_lw=1000   #km/s
em_uw=10000 #km/s  maximum width allowed for emision lines
emn_lw=10
emn_uw=400
abs_lw=0
abs_uw=800  #km/s  maximum width allowed for absortion lines
d_center=1000 #km/s red-blue shift allowed to the center of each line
d_center_n=200
em_la=0.0
em_ua=1.0
ab_la=-1.0
ab_ua=0.0


R_UV=3300  #
R_Vis=5400 #
R_IR=3890  #




#------------------SiIV-OIV complex---------------------------#
OIV=1402.34
OIV_a=0.05
OIV_uw= kms_to_wl(em_uw,OIV)
OIV_lw= kms_to_wl(em_lw,OIV)
OIV_w=OIV_uw
OIV_uc= OIV + kms_to_wl(d_center,OIV) 
OIV_lc= OIV - kms_to_wl(d_center,OIV)


SiIV= 1396.75
SiIV_a=0.05
SiIV_uw= kms_to_wl(em_uw,OIV)
SiIV_lw= kms_to_wl(em_lw,OIV)
SiIV_w=SiIV_uw
print SiIV_uw,SiIV_lw
SiIV_uc= SiIV + kms_to_wl(d_center,OIV) 
SiIV_lc= SiIV - kms_to_wl(d_center,OIV) 


delta=SiIV-OIV
#-------------------------------------------------------------#

guesses_SiIV = [OIV_a,OIV,OIV_w, SiIV_a, SiIV, SiIV_w,OIV_a,OIV,OIV_w/2, SiIV_a, SiIV, SiIV_w/2]
limits_SiIV=[(0,1), (OIV_lc , OIV_uc), (OIV_lw , OIV_uw), (0,1), (SiIV_lc , SiIV_uc), (SiIV_lw,SiIV_uw), (0,1), (OIV_lc , OIV_uc), (OIV_lw , OIV_uw/2), (0,1), (SiIV_lc , SiIV_uc), (SiIV_lw,SiIV_uw/2)]
limited_SiIV=[(True,False),(True,True), (True,True),   (True,False),  (False,False), (False,False),  
              (True,False),(False,False), (True,True), (True,False),  (False,False), (False,False)]
tied_SiIV=["","","",   "","","", 
           "","p[1]","",  "","p[4]",""]
xmin_SiIV=1360.0
xmax_SiIV=1440.0


#--------------------------- CIV complex -------------------------#
i=0
xmin_CIV=1460.0
xmax_CIV=1720.0



NIV=1486.50
NIV_a=0.1#sp.data( np.argmin( np.abs(NIV - sp.xarr)) ) 
NIV_uw= kms_to_wl(em_uw,NIV)/3.0
NIV_lw= kms_to_wl(em_lw,NIV)
NIV_uc= NIV + kms_to_wl(d_center,NIV)
NIV_lc= NIV - kms_to_wl(d_center,NIV)
NIV_w=NIV_uw


CIV=1550.7
CIV_a=2.0
CIV_uw=kms_to_wl(em_uw,CIV)
CIV_lw=kms_to_wl(em_lw,CIV)
CIV_uc=CIV + kms_to_wl(d_center,CIV)
CIV_lc=CIV - kms_to_wl(d_center,CIV)
CIV_w=CIV_uw

CIV1=1548.20
CIV1_a=2.0
CIV1_uw=kms_to_wl(em_uw,CIV)
CIV1_lw=kms_to_wl(em_lw,CIV)
CIV1_uc=CIV1 + kms_to_wl(d_center,CIV)
CIV1_lc=CIV1 - kms_to_wl(d_center,CIV)
CIV1_w=CIV1_uw




HeII=1640.72
HeII_a=0.5
HeII_uw=kms_to_wl(em_uw/3.0,HeII)
HeII_lw=kms_to_wl(em_lw,HeII)
HeII_uc=HeII + kms_to_wl(d_center,HeII)
HeII_lc=HeII - kms_to_wl(d_center,HeII)
HeII_w=HeII_uw


OIII= 1666.14
OIII_a=0.5
OIII_uw= kms_to_wl(em_uw/3.0,OIII)
OIII_lw= kms_to_wl(em_lw,OIII)
OIII_uc= OIII + kms_to_wl(d_center,OIII)
OIII_lc= OIII - kms_to_wl(d_center,OIII)
OIII_w=OIII_uw

OIII1=1660.80
OIII1_a=0.5
OIII1_uw= kms_to_wl(em_uw/3.0,OIII)
OIII1_lw= kms_to_wl(em_lw,OIII)
OIII1_uc= OIII1 + kms_to_wl(1.5*d_center,OIII)
OIII1_lc= OIII1 - kms_to_wl(1.5*d_center,OIII)
OIII1_w=OIII_uw

NIV=1718.55
NIV_a=0.5
NIV_uw= kms_to_wl(em_uw/2.0,NIV)
NIV_lw= kms_to_wl(em_lw,NIV)
NIV_uc= NIV + kms_to_wl(d_center,NIV)
NIV_lc= NIV - kms_to_wl(d_center,NIV)
NIV_w=NIV_uw

NIII=1754.00
NIII_a=0.5
NIII_uw= kms_to_wl(em_uw/3.0,NIII)
NIII_lw= kms_to_wl(em_lw,NIII)
NIII_uc= NIII + kms_to_wl(d_center,NIII)
NIII_lc= NIII - kms_to_wl(d_center,NIII)
NIII_w=NIII_uw


NIII1=1752.16
NIII1_a=0.5
NIII1_uw= kms_to_wl(em_uw/3.0,NIII)
NIII1_lw= kms_to_wl(em_lw,NIII)
NIII1_uc= NIII1 + kms_to_wl(1.5*d_center,NIII)
NIII1_lc= NIII1 - kms_to_wl(1.5*d_center,NIII)
NIII1_w=NIII1_uw


NIII2=1748.65
NIII2_a=0.5
NIII2_uw= kms_to_wl(em_uw/3.0,NIII)
NIII2_lw= kms_to_wl(em_lw,NIII)
NIII2_uc= NIII2 + kms_to_wl(1.5*d_center,NIII)
NIII2_lc= NIII2 - kms_to_wl(1.5*d_center,NIII)
NIII2_w=NIII2_uw
guessesNIII=[NIII_a,NIII,NIII_w]
limits_NIII=[(0,1), (NIII_lc,NIII_uc),(NIII_lw,NIII_uw)]
limited_NIII=[ (True,False),  (True,True), (True,True)]
tied_NIII=["","",""]
        

SiII=1818.17
SiII_a=0.5
SiII_uw= kms_to_wl(em_uw/3.0,SiII)
SiII_lw= kms_to_wl(em_lw,SiII)
SiII_uc= SiII + kms_to_wl(1.5*d_center,SiII)
SiII_lc= SiII - kms_to_wl(1.5*d_center,SiII)
SiII_w=SiII_uw
guessesSiII=[SiII_a,SiII,SiII_w]
limits_SiII=[(0,1), (SiII_lc,SiII_uc),(SiII_lw,SiII_uw)]
limited_SiII=[ (True,False),  (True,True), (True,True)]
tied_SiII=["","",""]


#-------------------------------------------------------------#

guesses_CIV=[NIV_a,NIV,NIV_w, CIV_a, CIV, CIV_w,CIV_a,CIV,CIV_w/3.0, CIV1_a, CIV1, CIV1_w, CIV1_a,CIV1,CIV1_w/3.0,   HeII_a,HeII,HeII_w,  HeII_a,HeII,HeII_w/3.0,   OIII_a,OIII,OIII_w, 0.4084*OIII1_a,OIII1,OIII1_w, NIII_a,NIII,NIII_w, 0.3415*NIII1_a,NIII1,NIII1_w, 1.09756*NIII2_a,NIII2,NIII2_w  ]
limits_CIV= [ (0,1), (NIV_lc , NIV_uc), (NIV_lw,NIV_uw), (0,1), (CIV_lc , CIV_uc), (CIV_lw,CIV_uw), (0,1), (CIV_lc , CIV_uc), (0,CIV_uw/3.0), (0,1), (CIV1_lc , CIV1_uc), (CIV1_lw,CIV1_uw), (0,1), (CIV1_lc , CIV1_uc), (0,CIV1_uw/2.0),   (0,1), (HeII_lc , HeII_uc), (HeII_lw/2.0,HeII_uw), (0,1), (HeII_lc , HeII_uc), (HeII_lw/2.0,HeII_uw/3.0),   (0,1), (OIII_lc,OIII_uc), (OIII_lw/2.0,OIII_uw/2.0),(0,1), (OIII1_lc,OIII1_uc), (OIII1_lw/2.0,OIII1_uw/2.0),(0,1), (NIII_lc,NIII_uc),(NIII_lw/2.0,NIII_uw/2.0),(0,1), (NIII1_lc,NIII1_uc),(NIII1_lw/2.0,NIII1_uw/2.0),(0,1), (NIII2_lc,NIII2_uc),(NIII2_lw/2.0,NIII2_uw/2.0), (0,1), (SiII_lc,SiII_uc),(SiII_lw,SiII_uw) ]
limited_CIV=[(True,False),(True,True), (True,True), (True,False),  (True,True), (True,True),   (True,False),(True,True), (True,True),  (True,False),  (True,True), (True,True),  (True,False),(True,True), (True,True),  (True,False),  (True,True), (True,True),  (True,False),(True,True), (True,True), (True,False),  (True,True), (True,True), (True,False),(True,True), (True,True), (True,False),(True,True), (True,True),  (True,False),(True,True), (True,True), (True,False),(True,True), (True,True),(True,False),(True,True), (True,True) ]    
tied_CIV=["", "", "",    "","","",  "","","",  "p[3]","p[4]+2.5","p[5]",    "p[6]","p[7]+2.5","p[8]",   "","","p[5]",  "","p[16]","",  "","","",    "0.4084*p[21]","p[22]-5.34","p[23]",   "", "", "",  "0.3415*p[27]", "p[28]-3.51", "p[29]",  "1.09756*p[27]", "p[28]-5.35", "p[29]",    "","",""]
xmin_CIV=1430
xmax_CIV=1800

#-----------------------CIII complex--------------------------#
xmin_CIII=1750.0
xmax_CIII=2100.0


CIII=1908.734
CIII_a=0.05
CIII_uw= kms_to_wl(em_uw/3.0,CIII)
CIII_lw= kms_to_wl(em_lw/2.0,CIII)
CIII_uc= CIII + kms_to_wl(d_center,CIII)
CIII_lc= CIII - kms_to_wl(d_center,CIII)
CIII_w=CIII_uw

SiIII=1892.03
SiIII_a=0.05
SiIII_uw= kms_to_wl(em_uw/3.0,CIII)
SiIII_lw= kms_to_wl(em_lw/2.0,CIII)
SiIII_uc= SiIII + kms_to_wl(d_center,CIII)
SiIII_lc= SiIII - kms_to_wl(d_center,CIII)
SiIII_w=SiIII_uw

AlIII1=1862.78
AlIII1_a=0.01
AlIII1_w=15.0
AlIII1_uw= kms_to_wl(em_uw/3.0,AlIII1)
AlIII1_lw= kms_to_wl(em_lw,AlIII1)
AlIII1_uc= AlIII1 + kms_to_wl(d_center,AlIII1)
AlIII1_lc= AlIII1 - kms_to_wl(d_center,AlIII1)
AlIII1_w=AlIII1_uw

AlIII2=1854.72
AlIII2_a=0.01
AlIII2_uw= kms_to_wl(em_uw/3.0,AlIII1)
AlIII2_lw= kms_to_wl(em_lw,AlIII1)
AlIII2_uc= AlIII2 + kms_to_wl(d_center,AlIII1)
AlIII2_lc= AlIII2 - kms_to_wl(d_center,AlIII1)
AlIII2_w=AlIII2_uw





#-------------------------------------------------------------#

guesses_CIII=[CIII_a,CIII,CIII_w,CIII_a,CIII,CIII_w/3.0,SiIII_a,SiIII,SiIII_w,SiIII_a,SiIII,SiIII_w/3.0, AlIII1_a, AlIII1,AlIII1_w,AlIII2_a, AlIII2,AlIII2_w] 
limits_CIII=[ (0,1), (CIII_lc , CIII_uc), (CIII_lw,CIII_uw),(0,1), (CIII_lc , CIII_uc), (CIII_lw,CIII_uw/3.0), (0,1), (SiIII_lc , SiIII_uc), (SiIII_lw,SiIII_uw), (0,1), (SiIII_lc , SiIII_uc), (SiIII_lw,SiIII_uw/3.0), (0,1), (AlIII1_lc,AlIII1_uc), (AlIII1_lw,AlIII1_uw), (0,1), (AlIII2_lc,AlIII2_uc), (AlIII2_lw,AlIII2_uw)] 
limited_CIII=[ (True,False),(True,True), (True,True),     (True,False),  (True,True), (True,True),  (True,False),  (True,True), (True,True), (True,False),  (True,True), (True,True),(True,False),  (True,True), (True,True), (True,False),  (True,True), (True,True)]
tied_CIII=[ "","","",    "", "p[1]","",   "p[0]","p[1] - 16.704","p[2]",   "","p[7]","p[5]",   "","p[1] - 54.014","",   "p[12]","p[13]+8.06","p[14]"]


guesses_C=[NIV_a,NIV,NIV_w, CIV_a, CIV, CIV_w,
           CIV_a,CIV,CIV_w/2.0, CIV1_a, CIV1, CIV1_w,
           CIV1_a,CIV1,CIV1_w/2.0,   HeII_a,HeII,HeII_w/2.0,
           HeII_a,HeII,HeII_w/3.0,   OIII_a,OIII,OIII_w/2.0,
           0.4084*OIII1_a,OIII1,OIII1_w, NIII_a,NIII,NIII_w/2.0,
           0.3415*NIII1_a,NIII1,NIII1_w/2.0, 1.09756*NIII2_a,NIII2,NIII2_w/2.0,
           CIII_a,CIII,CIII_w,CIII_a,CIII,CIII_w/3.0,
           SiIII_a,SiIII,SiIII_w,SiIII_a,SiIII,SiIII_w/3.0, 
           AlIII1_a, AlIII1,AlIII1_w,AlIII2_a, AlIII2,AlIII2_w,
           SiII_a, SiII,SiII_w, NIV_a,NIV,NIV_w/2.0] 



limits_C= [ (0,1), (NIV_lc , NIV_uc), (NIV_lw,NIV_uw), (0,1), (CIV_lc , CIV_uc), (CIV_lw,CIV_uw), 
            (0,1), (CIV_lc , CIV_uc), (0,CIV_uw/2.0), (0,1), (CIV1_lc , CIV1_uc), (CIV1_lw,CIV1_uw),
            (0,1), (CIV1_lc , CIV1_uc), (0,CIV1_uw/2.0),   (0,1), (HeII_lc , HeII_uc), (HeII_lw,HeII_uw),
            (0,1), (HeII_lc , HeII_uc), (HeII_lw,HeII_uw/3.0),   (0,1), (OIII_lc,OIII_uc), (OIII_lw,OIII_uw),
            (0,1), (OIII1_lc,OIII1_uc), (OIII1_lw,OIII1_uw),(0,1), (NIII_lc,NIII_uc),(NIII_lw,NIII_uw),
            (0,1),(NIII1_lc,NIII1_uc),(NIII1_lw,NIII1_uw),(0,1), (NIII2_lc,NIII2_uc),(NIII2_lw,NIII2_uw), 
            (0,1), (CIII_lc , CIII_uc), (CIII_lw,CIII_uw),(0,1), (CIII_lc , CIII_uc), (CIII_lw,CIII_uw/3.0),
            (0,1), (SiIII_lc , SiIII_uc), (SiIII_lw,SiIII_uw), (0,1), (SiIII_lc , SiIII_uc), (SiIII_lw,SiIII_uw/3.0), 
            (0,1), (AlIII1_lc,AlIII1_uc), (AlIII1_lw,AlIII1_uw), (0,1), (AlIII2_lc,AlIII2_uc), (AlIII2_lw,AlIII2_uw),
            (0,1), (SiII_lc,SiII_uc),(SiII_lw,SiII_uw), (0,1), (NIV_lc , NIV_uc), (0,NIV_uw/2.0)] 


limited_C=[(True,False),(True,True), (True,True),        (True,False),  (True,True), (True,True),    
           (True,False),(True,True), (True,True),        (True,False),(True,True), (True,True),
           (True,False),(True, True),(True,True),        (True,False),(True,True), (True,True),   
           (True,False),  (True,True), (True,True),      (True,False),(True,True), (True,True),    
           (True,False),  (True,True), (True,True),      (True,False),(True,True), (True,True), 
           (True,False),(True,True), (True,True),        (True,False),(True,True), (True,True), 
           (True,False),(True,True), (True,True),        (True,False),(True,True), (True,True),    
           (True,False),  (True,True), (True,True),      (True,False),  (True,True), (True,True),
           (True,False),  (True,True), (True,True),      (True,False),  (True,True), (True,True),
           (True,False),  (True,True), (True,True),      (True,False),  (True,True), (True,True)]    
tied_C=["", "", "",                "","","",
        "", "", "",                "","p[4]-2.5","p[5]", 
        "","p[7]-2.5","p[8]",   "","","", 
        "","p[16]","",              "","","",  
        "0.4084*p[21]","p[22]-5.34","p[23]",   "", "", "",  
        "0.3415*p[27]", "p[28]-3.51", "p[29]",  "1.09756*p[27]", "p[28]-5.35", "p[29]", 
        "","","",    "", "p[37]","", 
        "","","p[38]",   "","p[43]","",  
        "","","",   "p[48]","p[49]-8.06","p[50]",
        "", "", "", "", "", "" ]


xmin_C=1430.0
xmax_C=2000.0



#------------------MgII1-OIV complex---------------------------#

MgII2=2802.71
MgII2_a=0.5
MgII2_uw= kms_to_wl(em_uw/2,MgII2)
MgII2_lw= kms_to_wl(em_lw,MgII2)
MgII2_uc= MgII2 + kms_to_wl(d_center,MgII2)
MgII2_lc= MgII2 - kms_to_wl(d_center,MgII2)
MgII2_w=MgII2_uw



MgII1= 2795.53
MgII1_a=0.5
MgII1_uw= kms_to_wl(em_uw/2,MgII2)
MgII1_lw= kms_to_wl(em_lw,MgII2)
MgII1_uc= MgII1 + kms_to_wl(d_center,MgII2)
MgII1_lc= MgII1 - kms_to_wl(d_center,MgII2)
MgII1_w=MgII1_uw
delta=MgII1-MgII2
#-------------------------------------------------------------#

guesses_MgII = [MgII2_a,MgII2,MgII2_w,MgII2_a,MgII2,MgII2_w/2.0, MgII1_a, MgII1, MgII1_w, MgII1_a, MgII1, MgII1_w/2.0 ]
limits_MgII=[(0,1), (MgII2_lc , MgII2_uc), (MgII2_lw , MgII2_uw), (0,1), (MgII2_lc , MgII2_uc), (MgII2_lw , MgII2_uw), (0,1), (MgII1_lc , MgII1_uc), (MgII1_lw,MgII1_uw) , (0,1), (MgII1_lc , MgII1_uc), (MgII1_lw,MgII1_uw) ]
limited_MgII=[(True,False),(True,True), (True,True), (True,False),  (True,True), (True,True),(True,False),(True,True), (True,True), (True,False),  (True,True), (True,True) ]
tied_MgII=["","","", "","","",  "","p[1]-7.18","p[2]",  "","p[1]-7.18","p[5]"]
xmin_MgII=2700
xmax_MgII=2900
#7.18



#--------------------------- HeI complex -------------------------#


HeI=5875.50
HeI_a=0.02#sp.data( np.argmin( np.abs(HeI - sp.xarr)) ) 
HeI_uw= kms_to_wl(em_uw,HeI)
HeI_lw= kms_to_wl(em_lw,HeI)
HeI_uc= HeI + kms_to_wl(d_center,HeI)
HeI_lc= HeI - kms_to_wl(d_center,HeI)
HeI_w=HeI_uw

#----------------------------------------------------------------#

guesses_HeI = [ HeI_a, HeI, HeI_w, HeI_a, HeI, HeI_w/3]
limits_HeI=[(0,1), (HeI_lc , HeI_uc), (HeI_lw,HeI_uw),  (0,1), (HeI_lc , HeI_uc), (HeI_lw,HeI_uw/2)]
limited_HeI=[(True,False),(True,True), (True,True), (True,False),  (True,True), (True,True)]
tied_HeI=["","","",  "","p[1]",""]
xmin_HeI=5820.5
xmax_HeI=5920.0


#--------------------------- Hbeta complex -------------------------#


HeIIb=4685.65
HeIIb_a=0.05
HeIIb_uw=kms_to_wl(em_uw,HeIIb)
HeIIb_lw=kms_to_wl(em_lw,HeIIb)
HeIIb_uc=HeIIb + kms_to_wl(d_center,HeIIb)
HeIIb_lc=HeIIb - kms_to_wl(d_center,HeIIb)
HeIIb_w=HeIIb_uw


Hbeta=4861.32
Hbeta_a=0.1
Hbeta_uw= kms_to_wl(em_uw,Hbeta)
Hbeta_lw= kms_to_wl(em_lw,Hbeta)
Hbeta_uc= Hbeta + kms_to_wl(d_center,Hbeta)
Hbeta_lc= Hbeta - kms_to_wl(d_center,Hbeta)
Hbeta_w=Hbeta_uw

Hbetan=4861.32
Hbetan_a=0.1
Hbetan_uw=kms_to_wl(emn_uw,Hbetan)
Hbetan_lw=kms_to_wl(emn_lw,Hbetan)
Hbetan_uc=Hbetan + kms_to_wl(d_center,Hbetan)
Hbetan_lc=Hbetan - kms_to_wl(d_center,Hbetan)
Hbetan_w=Hbetan_uw

OIIIb=5006.84
OIIIb_a=0.03
OIIIb_uw= kms_to_wl(emn_uw,OIIIb)
OIIIb_lw= kms_to_wl(emn_lw,OIIIb)
OIIIb_uc= OIIIb + kms_to_wl(d_center,OIIIb)
OIIIb_lc= OIIIb - kms_to_wl(d_center,OIIIb)
OIIIb_w=OIIIb_uw

OIIIb1=4958.91
OIIIb1_a=0.03
OIIIb1_uw= kms_to_wl(emn_uw,OIIIb)
OIIIb1_lw= kms_to_wl(emn_lw,OIIIb)
OIIIb1_uc= OIIIb1 + kms_to_wl(d_center,OIIIb)
OIIIb1_lc= OIIIb1 - kms_to_wl(d_center,OIIIb)
OIIIb1_w=OIIIb1_uw
#-------------------------------------------------------------#

guesses_Hbeta=[Hbeta_a,Hbeta,Hbeta_w,    Hbeta_a,Hbeta,Hbeta_w/3.0,
               HeIIb_a,HeIIb,HeIIb_w,    OIIIb_a,OIIIb,OIIIb_w,
               OIIIb1_a,OIIIb1,OIIIb1_w, Hbetan_a,Hbetan,Hbetan_w]
limits_Hbeta= [ (0,1), (Hbeta_lc , Hbeta_uc), (Hbeta_lw,Hbeta_uw),   (0,1), (Hbeta_lc , Hbeta_uc), (Hbeta_lw,Hbeta_uw/3.0), 
                (0,1), (HeIIb_lc , HeIIb_uc), (HeIIb_lw,HeIIb_uw),   (0,1), (OIIIb_lc,OIIIb_uc), (OIIIb_lw,OIIIb_uw), 
                (0,1), (OIIIb1_lc,OIIIb1_uc), (OIIIb1_lw,OIIIb1_uw), (0,1), (Hbetan_lc , Hbetan_uc), (Hbetan_lw,Hbetan_uw)]
limited_Hbeta=[(True,False),(True,True), (True,True), (True,False),  (True,True), (True,True),
               (True,False),(True,True), (True,True), (True,False),  (True,True), (True,True),
               (True,False),(True,True), (True,True),(True,False),(True,True), (True,True)]    
tied_Hbeta=["", "", "",    "","","", 
            "","","",  "","","",   
            "0.33*p[9]","p[10]- 47.93","p[11]",   "","",""]
xmin_Hbeta=4620
xmax_Hbeta=5120



#--------------------------- H-ALPHA complex -------------------------#
em_lw=1000   #km/s
em_uw=10000 #km/s  maximum width allowed for emision lines
emn_lw=10
emn_uw=400
abs_lw=0
abs_uw=800  #km/s  maximum width allowed for absortion lines
d_center=4000 #km/s red-blue shift allowed to the center of each line
d_center_n=200


NII=6583.39
NII_a=0.05
NII_uw=kms_to_wl(emn_uw,NII)
NII_lw=kms_to_wl(emn_lw,NII)
NII_uc=NII + kms_to_wl(d_center,NII)
NII_lc=NII - kms_to_wl(d_center,NII)
NII_w=NII_uw


NII1=6548.06
NII1_a=0.05
NII1_uw=kms_to_wl(emn_uw,NII)
NII1_lw=kms_to_wl(emn_lw,NII)
NII1_uc=NII1 + kms_to_wl(d_center,NII)
NII1_lc=NII1 - kms_to_wl(d_center,NII)
NII1_w=NII1_uw



Halpha=6562.80
Halpha_a=0.1
Halpha_uw= kms_to_wl(em_uw,Halpha)
Halpha_lw= kms_to_wl(em_lw,Halpha)
Halpha_uc= Halpha + kms_to_wl(d_center,Halpha)
Halpha_lc= Halpha - kms_to_wl(d_center,Halpha)
Halpha_w=Halpha_uw

Halphan=6562.80
Halphan_a=0.1
Halphan_uw=kms_to_wl(emn_uw,Halphan)
Halphan_lw=kms_to_wl(emn_lw,Halphan)
Halphan_uc=Halphan + kms_to_wl(d_center,Halphan)
Halphan_lc=Halphan - kms_to_wl(d_center,Halphan)
Halphan_w=Halphan_uw

SII=6730.85
SII_a=0.03
SII_uw= kms_to_wl(emn_uw,SII)
SII_lw= kms_to_wl(emn_lw,SII)
SII_uc= SII + kms_to_wl(d_center,SII)
SII_lc= SII - kms_to_wl(d_center,SII)
SII_w=SII_uw

SII1=6716.47
SII1_a=0.03
SII1_uw= kms_to_wl(emn_uw,SII)
SII1_lw= kms_to_wl(emn_lw,SII)
SII1_uc= SII1 + kms_to_wl(d_center,SII)
SII1_lc= SII1 - kms_to_wl(d_center,SII)
SII1_w=SII1_uw

#-------------------------------------------------------------#

guesses_Halpha=[Halpha_a,Halpha,Halpha_w, Halpha_a,Halpha,Halpha_w/3.0,
                NII_a,NII,NII_w,          NII1_a,NII1,NII1_w,
                SII_a,SII,SII_w,         SII1_a,SII1,SII1_w,
                Halphan_a,Halphan,Halphan_w]
limits_Halpha= [ (0,1), (Halpha_lc , Halpha_uc), (Halpha_lw,Halpha_uw),  (0,1), (Halpha_lc , Halpha_uc), (Halpha_lw,Halpha_uw/3.0),
                 (0,1), (NII_lc , NII_uc), (NII_lw,NII_uw),              (0,1), (NII1_lc , NII1_uc), (NII1_lw,NII1_uw),
                 (0,1), (SII_lc,SII_uc), (SII_lw,SII_uw),                (0,1), (SII1_lc,SII1_uc), (SII1_lw,SII1_uw),
                 (0,1),(Halphan_lc , Halphan_uc), (Halphan_lw,Halphan_uw)]
limited_Halpha=[(True,False),(True,True), (True,True), (True,False),  (True,True), (True,True),
                (True,False),(True,True), (True,True), (True,False),  (True,True), (True,True),
                (True,False),(True,True), (True,True),  (True,False),(True,True), (True,True)  ,
                (True,False),(True,True), (True,True)] 
tied_Halpha=["", "", "",    "","","",
             "","","", "3*p[6]","p[7]-35.33","p[8]", 
             "","","",    "p[12]","p[13]-14.38" ,"",
             "","",""]
xmin_Halpha=6350
xmax_Halpha=6800


