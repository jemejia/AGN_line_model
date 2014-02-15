
from matplotlib import pylab
import numpy as np
import pyspeckit
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

#from pyspeckit.spectrum.models.template import template_fitter as temp_fitter
xmin=1250
xmax=9000
n_spec=1

def continuous_substraction( spec_number, magorder=-16):
    z=np.loadtxt("redshift_list.txt")
    model=np.loadtxt("./model/cont_model_" + str(spec_number) + ".txt")

    Lmodel=fdata=np.interp(sp.xarr,model[:,0],model[:,1])
    print "model=", model[:,1]
    print "Lmodel=", Lmodel
    print "wl_model=", model[:,0]
    print "xarr", sp.xarr
    z=z[spec_number]
    dL= cosmo.luminosity_distance(z)*3.08567758e24 
    print "dL",dL
    dl2=dL*dL
    print "dl2",dl2
    Fmodel=np.divide( Lmodel, 4*np.pi*dl2 )
    magorder=-1.0*magorder + 9
    Fmodel=4*Fmodel*np.power(10, magorder ) -0.15
    Fmodel=4*Fmodel*np.power(10, magorder ) - 0.025 
    Fmodel=np.multiply(Lmodel, np.power(10,magorder) )
    print Fmodel
    print sp.data

    
    pylab.figure()
    pylab.ylabel(r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.plot(sp.xarr,sp.data,'k')
    #pylab.plot(sp.xarr,Fmodel,'r')
    pylab.show()
    Lmodel=4*np.pi*Fmodel* dl2
    sp.data=sp.data - Fmodel 
    return Fmodel,Lmodel

def smooth(data,smooth,smoothtype='gaussian',downsample=False,downsample_factor=None,
        convmode='same'):
    """
    Smooth and downsample the data array

    Parameters
    ----------
    smooth  :  float 
        Number of pixels to smooth by
    smoothtype : [ 'gaussian','hanning', or 'boxcar' ]
        type of smoothing kernel to use
    downsample :  bool 
        Downsample the data?
    downsample_factor  :  int 
        Downsample by the smoothing factor, or something else?
    convmode : [ 'full','valid','same' ]
        see :mod:`numpy.convolve`.  'same' returns an array of the same length as
        'data' (assuming data is larger than the kernel)
    """
    
    roundsmooth = smooth # can only downsample by integers

    if downsample_factor is None and downsample:
        downsample_factor = roundsmooth
    elif downsample_factor is None:
        downsample_factor = 1
    downsample_factor = 1
    if smooth > len(data) or downsample_factor > len(data):
        raise ValueError("Error: trying to smooth by more than the spectral length.")

    
    
    xkern  = np.linspace(-5*smooth,5*smooth,smooth*11)
    kernel = np.exp(-xkern**2/(2*(smooth/np.sqrt(8*np.log(2)))**2))
    kernel /= kernel.sum()
    if len(kernel) > len(data):
        lengthdiff = len(kernel)-len(data)
        if lengthdiff % 2 == 0: # make kernel same size as data
            kernel = kernel[lengthdiff/2:-lengthdiff/2]
        else: # make kernel 1 pixel smaller than data but still symmetric
            kernel = kernel[lengthdiff/2+1:-lengthdiff/2-1]
    
    # deal with NANs or masked values
    if hasattr(data,'mask'):
        if type(data.mask) is np.ndarray:
            OK = True - data.mask
            if OK.sum() > 0:
                data = pyspeckit.interpolation._interp(np.arange(len(data)),np.arange(len(data))[OK],data[OK])
            else:
                data = OK
    if np.any(True - np.isfinite(data)):
        OK = np.isfinite(data)
        if OK.sum() > 0:
            data = pyspeckit.interpolation._interp(np.arange(len(data)),np.arange(len(data))[OK],data[OK])
        else:
            data = OK

    if np.any(True - np.isfinite(data)):
        raise ValueError("NANs in data array after they have been forcibly removed.")

    smdata = np.convolve(data,kernel,convmode)#[::downsample_factor]

    return smdata







def kms_to_wl(kms,line_center):
    c=3e5
    wl=kms*line_center/c
    return wl


    
def blue_bump_fitter(line_name,spec_number,xmin,xmax,fe_template,balmer_template):
    
    #-----------------------------defining fe template fitter ------------------------------------#
    #template_file="./templates/fe2_UV.dat"
    #template = pyspeckit.Spectrum(template_file)

    for i in range(5):
        
        def model(xarr, scale_fe,shift_fe, smooth_factor_fe):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
        

            new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                       fe_template.xarr.as_unit('angstroms')+shift_fe, 
                                                                       fe_template.data )
        
            #data=smooth(new_fe_template.data,smooth_factor_fe , downsample=False    )
            #new_fe_template.data=data
            #new_balmer_template = scale_balmer*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
            #                                         balmer_template.xarr.as_unit('angstroms'), 
            #                                         balmer_template.data)
            

            new_template=  new_fe_template #+ new_balmer_template 
            return new_template


        


        modelclass = pyspeckit.spectrum.models.model.SpectralModel(model,
                                                                   3, 
                                                                   parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                   parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                   parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                   shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                   centroid_par='shift_fe' )
        modelclass.__name__ = "template"



        #-----------------------------defining fe template fitter ------------------------------------#

        
        #-----------restarting plotter------------#
        sp.plotter.refresh()    
        sp.plotter.reset_limits()
        sp.plotter(xmin=xmin,xmax=xmax)
        #-----------restarting plotter------------#
        
        
        # -----------------  fitting ---------------------
        exclude_file = "./excludes/exclude_" + line_name + ".txt"
        if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
            exclude_data=np.loadtxt(exclude_file,skiprows=2)
            exclude=exclude_data[spec_number,:]
        else:
            exclude=[]
            
            
        
        sp.Registry.add_fitter('template',modelclass,3,multisingle='multi')
        sp.specfit(fittype='template',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax)
        print "sp1 lenght=", len(sp.data)
        sp.specfit.fullsizemodel()  
        print "sp2 lenght=", len(sp.data)
        fe_template.data*=sp.specfit.parinfo.values[0]
        fe_template.xarr+=sp.specfit.parinfo.values[1]
        new_template=fe_template.copy()
        
        def suavizado(xarr, smooth_factor_fe,template=new_template.copy()):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
            
            
            data=smooth(template.data,smooth_factor_fe , downsample=False)   
            template.data=data
            n_template =pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                        template.xarr.as_unit('angstroms'), 
                                                        template.data)
            fe_template.data=data
            #new_template=  fe_template #+ new_balmer_template 
            
            
            return n_template
    
        
        soft = pyspeckit.spectrum.models.model.SpectralModel(suavizado,
                                                             1, 
                                                             parnames=['smooth_fe'],
                                                             parlimited=[(True,False)],
                                                             parlimits=[(0,0)],
                                                             shortvarnames=(r'sigma_x') )
        soft.__name__ = "suavizado"
        
        sp.Registry.add_fitter('suavizado',soft,1,multisingle='multi')
        sp.specfit(fittype='suavizado',exclude=exclude,guesses=[1.0],xmin=xmin,xmax=xmax)
 
        # -----------------  fitting ---------------------
            
    # -----------------  ploting  fit---------------------

    
    fitmin=sp.specfit.xmin
    fitmax=sp.specfit.xmax
    sp.specfit.fullsizemodel() 
    plot_file="./plots/fe_fit_" + line_name + "_" + str(spec_number) + ".png"
    pylab.figure()
    pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=2*sp.specfit.residuals[fitmin:fitmax].max())
    pylab.xlim(xmin=xmin-300,xmax=xmax)
    pylab.ylabel(r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
    pylab.savefig(plot_file)
    # -----------------  ploting fit---------------------
    
    
    # ----------plotting residuals-------------
    
    pylab.figure()
    pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals[fitmin:fitmax].max())
    pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
    pylab.ylabel(r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    print len(sp.xarr[fitmin:fitmax]), len(sp.specfit.residuals[fitmin:fitmax])
    pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
    plot_file="./plots/" + line_name +"_res_" + str(spec_number) + ".png"
    
    pylab.savefig(plot_file)
    # ----------plotting residuals-------------#
    
    #----------redifining spectrum ------------------#
    sp.data=sp.data - sp.specfit.fullmodel
    #----------redifining spectrum ------------------#
    return sp.specfit.fullmodel
        
def line_fitter(line_name, spec_number,guesses, limits, limited, tied, xmin, xmax, offset=-0.4, magorder=-16):
    #-----------restarting plotter------------#
    sp.plotter.refresh()    
    sp.plotter.reset_limits()
    sp.plotter(xmin=xmin,xmax=xmax)
    #-----------restarting plotter------------#

    # -----------------  fitting ---------------------
    exclude_file = "./excludes/exclude_" + line_name + ".txt"
    if os.path.exists(exclude_file) and os.path.isfile(exclude_file):
        exclude_data=np.loadtxt(exclude_file,skiprows=2)
        exclude=exclude_data[spec_number,:]
    else:
        exclude=[]
    
    sp.specfit(exclude=exclude,guesses = guesses, limits=limits, limited=limited, tied=tied, annotate = False,fittype='gaussian',xmin=xmin,xmax=xmax,multifit=False)
    # -----------------  fitting ---------------------
   
    # ----------------- fit plotting  ---------------------
    sp.specfit.plot_components(add_baseline=False)
    
    sp.plotter.refresh()
    sp.specfit.plot_fit(annotate=False)
    # -----------------  fit plotting  ---------------------

    # -----------------  line measuring  ---------------------
    sp.measure(z = 1.55, fluxnorm = 10**(-1.0*mag_order))
    # ----------------- TO DEVELOP!!!!  ---------------------
    # ----------------- SiOIV line measuring  ---------------------

    # ----------plotting residuals-------------
    plot_file="./plots/" + line_name + "_" + str(spec_number) + ".png"
    sp.plotter.figure.savefig(plot_file)
    
    fitmin=sp.specfit.xmin
    fitmax=sp.specfit.xmax
    pylab.figure()
    pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals.max())
    pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
    pylab.ylabel(r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals,'k')
    plot_file="./plots/" + line_name +"_res_" + str(spec_number) + ".png"
    
    pylab.savefig(plot_file)
    # ----------plotting residuals-------------#
    
    
    #---------------measuring----------------------#
    sp.measure(z = 1.51831, fluxnorm = 10**magorder)

    for i, line in enumerate(sp.measurements.lines.keys()):

        # If this line is not in our database of lines, don't try to annotate it
        if line not in sp.speclines.optical.lines.keys(): continue

        x = sp.measurements.lines[line]['modelpars'][1]   # Location of the emission line


    # Print out spectral line information
    filename="./plots/" + line_name + "_" + str(spec_number) + ".txt"
    f=open(filename,"w")

    
    for line in sp.measurements.lines.keys():
        f.write("#Line   Flux (erg/s/cm^2)     Amplitude (erg/s/cm^2)    FWHM (Angstrom)   Luminosity (erg/s)")
        print line, sp.measurements.lines[line]['flux'], sp.measurements.lines[line]['amp'], sp.measurements.lines[line]['fwhm'], sp.measurements.lines[line]['lum']
        print_line = '{0}  {1}  {2}  {3} {4} {5}\n'.format(line, sp.measurements.lines[line]['flux'],sp.measurements.lines[line]['amp'],sp.measurements.lines[line]['fwhm'],sp.measurements.lines[line]['lum'], sp.measurements.lines[line]['modelpars'][1] )
        f.write(print_line)
    f.close()
    return sp.specfit.fullmodel
    #---------------measuring----------------------#
pylab.ion()



em_lw=1000   #km/s
em_uw=10000 #km/s  maximum width allowed for emision lines
emn_lw=10
emn_uw=400
abs_lw=0
abs_uw=800  #km/s  maximum width allowed for absortion lines
d_center=4000 #km/s red-blue shift allowed to the center of each line
d_center_n=00
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
xmin_SiIV=1350.0
xmax_SiIV=1450.0


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






for i in range(n_spec,13):
    

    #------------------- -------Shared features for fitting. Finished -----------------------#
    #------------opening spectrum i---------------------#
    spectrum_file='../AGN_spectra/spectrum' + str(i) + '.txt'
    sp = pyspeckit.Spectrum(spectrum_file)
    #sp.smooth(2)
    #------------opening spectrum i---------------------#
    
    
    # -----------set up units properly------------------#
    mag_order=np.int((-1)*np.round(np.log10(np.mean(sp.data))))
    sp.xarr.units='angstroms'
    sp.xarr.xtype = 'wavelength'
    sp.units = r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$'
    sp.data *= 10**(mag_order)    
    #-------------- set up units properly------------#
    copy=sp.copy()
    copy2=sp.copy()
    copy1=sp.copy()
    copy3=copy.copy()
    #----------starting plotter--------------#
    copy.plotter(xmin=xmin,xmax=6900,ymin=-1.3*np.abs(copy2.data.min()),ymax=1.1*copy1.data.max())
    copy3.plotter(xmin=xmin,xmax=6900,ymin=-1.3*np.abs(copy2.data.min()),ymax=1.1*copy1.data.max())
    #----------starting  plotter--------------#


    #-------------continuous fitting-------------------#
    exclude_file = "./excludes/exclude_cont.txt"
    exclude_cont=np.loadtxt(exclude_file,skiprows=2)

    limit_file = "./cont_limits.txt"
    cont_limits=np.loadtxt(limit_file)
    slope_break=4230#cont_limits[:,1]
    slope_lim=4200#cont_limits[:,0]
    wlmin=sp.xarr[0]
    wlmax=sp.xarr[-1]
    
    
    pylab.xscale('log')#,subsx=[3,4])
    """
    copy.baseline.powerlaw=True
    copy.baseline(xmin=2100, xmax=xmax, exclude=exclude_cont[i,:], subtract=False, reset_selection=False, highlight_fitregion=True,powerlaw=True,quiet=False,LoudDebug=False,annotate=True)
    
    copy.crop(2200.0,copy.xarr[-1])
    copy3.baseline.powerlaw=True
    copy3.baseline(xmin=xmin, xmax=2250, exclude=exclude_cont[i,:], subtract=False, reset_selection=False, highlight_fitregion=True,powerlaw=True,quiet=False,LoudDebug=False,annotate=True)
    copy3.crop(copy3.xarr[0],2199.6)
    
    copy1.crop(1400,1600)
    pylab.ylim(ymin=-1.3*np.abs(copy1.data.min()),ymax=1.1*copy1.data.max())
    pylab.xlim(xmin=xmin,xmax=6900)
    pylab.ylabel(r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.xscale('log')
    pylab.plot(copy.xarr,copy.data,'k')
    pylab.plot(copy3.xarr,copy3.data,'k')
    pylab.plot(copy.xarr,copy.baseline.basespec,'b')
    pylab.plot(copy3.xarr,copy3.baseline.basespec,'r')
    plot_file="./plots/break_continuous_" + str(i) + ".png"
    copy.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)
    
    
    break_pos=np.argmin( np.abs( sp.xarr - slope_break) )

    copy_data=copy.data - copy.baseline.basespec
    copy3_data=copy3.data - copy3.baseline.basespec
    data=np.concatenate((copy_data,copy3_data))
    xarr=np.concatenate((copy.xarr,copy3.xarr))
    pylab.show()
    #dat=sp.data - copy.baseline.basespec
    sp=pyspeckit.Spectrum(xarr=xarr,data=data)
    sp.xarr.units='angstroms'
    sp.xarr.xtype = 'wavelength'
    sp.units = r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$'
    #sp.plotter.refresh()    
    sp.plotter(xmin=xmin,xmax=xmax,ymin=1.1*copy2.data.min(),ymax=3.5*copy2.data.max())
    
    
    plot_file="./plots/continuous_" + str(i) + ".png"
    copy.plotter.figure.savefig(plot_file)
    """
    #pylab.xscale('log',subsx=[3,4])
    
    sp.plotter(xmin=xmin,xmax=xmax,ymin=1.1*copy2.data.min(),ymax=3.5*copy2.data.max())
    #sp.baseline.powerlaw=True
    #sp.baseline(xmin=xmin, xmax=xmax, exclude=exclude_cont[i,:], subtract=True, reset_selection=False, highlight_fitregion=False,powerlaw=True,quiet=False,LoudDebug=False,annotate=False)

    continuous,L_cont=continuous_substraction( i, magorder=-1*mag_order)
    copy1.crop(1400,1600)
    pylab.ylim(ymin=-1.3*np.abs(copy1.data.min()),ymax=1.1*copy1.data.max())
    pylab.xlim(xmin=xmin,xmax=6900)
    pylab.ylabel(r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.xscale('log')
    pylab.plot(sp.xarr,sp.data,'k')
    plot_file="./plots/no_continuous_" + str(i) + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)
    
    #-------------continuous fitting-------------------#

    
    #------------------- -------Shared features for fitting. Finished ---------------------------#
    


    
    
    
    xmin_CIII=1940
    xmax_CIII=2010

    #fe_fitter("C",i,xmin_CIII,xmax_CIII)
    #line_fitter("C", i,guesses_C, limits_C, limited_C, tied_C, xmin_C, xmax_C,magorder=-1*mag_order)

    #line_fitter("HeI", i, guesses_HeI, limits_HeI, limited_HeI, tied_HeI, xmin_HeI, xmax_HeI)
    
    
    xmin_bb=3640
    xmax_bb=3650
    lambda_0=3675
    index=np.argmin( np.abs( lambda_0 - sp.xarr) )
    flux=np.mean(sp.data[index -30: index + 30])
    
    #blue_bump_fitter(xmin_bb,xmax_bb)
    #sp.data=sp.data - sp.specfit.fullmodel
    template_file="./templates/fe2_UV.dat"
    fe_template = pyspeckit.Spectrum(template_file)
    fe_template.xarr.units='angstroms'
    fe_template.xarr.xtype='wavelength'
    fe_template.crop(1555.0,fe_template.xarr[-1],units="angstroms")

    lambda_fe=2700
    index_sp=np.argmin( np.abs( lambda_fe - sp.xarr) )
    index_fe=np.argmin( np.abs( lambda_fe - fe_template.xarr) )
    flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
    fe_template.data*=flux/fe_template.data[index_fe]
    #fe_template.data*=np.mean(sp.data)/np.mean(fe_template.data)

    template_file="./templates/balmer_template2.txt"
    balmer_template = np.loadtxt(template_file)
    balmer_data=np.interp(sp.xarr,balmer_template[:,0], balmer_template[:,1])
    balmer_template=pyspeckit.Spectrum(data=balmer_data,xarr=sp.xarr)
    balmer_template.xarr.units='angstroms'
    balmer_template.xarr.xtype='wavelength'
    lambda_balmer=3645
    index=np.argmin( np.abs( lambda_balmer - sp.xarr) )
    flux=np.mean(sp.data[ index - 30 : index + 30 ])
    balmer_template.data*=flux/balmer_template.data[-1]
    balmer_data=np.zeros(len(balmer_template.data))
    balmer_data[:index]=balmer_template.data[:index]
    xmin_bb=3640
    xmax_bb=3650
    
    
    backup=sp.copy()
    sp.data[:index]=sp.data[:index] - balmer_template.data[:index]
    
    xmin_fe=2200
    xmax_fe=3080
    pylab.figure()
    pylab.plot(backup.xarr,backup.data)
    pylab.plot(balmer_template.xarr[:index],balmer_template.data[:index])
    pylab.plot(sp.xarr,sp.data)
    pylab.savefig("balmer_fit_"+str(i)+".pdf")
    backup=sp.copy()
    fe_fit=blue_bump_fitter("SBB",i,xmin_fe,xmax_fe,fe_template, balmer_template)
    Mg_fit=line_fitter("MgII", i, guesses_MgII, limits_MgII, limited_MgII, tied_MgII, xmin_MgII, xmax_MgII,offset=-0.05, magorder=-1*mag_order)
    print "FWHM_Mg=", sp.specfit.measure_approximate_fwhm()
    #fe_fitter("C",i,xmin_CIII,xmax_CIII)
    Si_fit=line_fitter("SiIV", i,guesses_SiIV, limits_SiIV, limited_SiIV, tied_SiIV, xmin_SiIV, xmax_SiIV, magorder=-1*mag_order)
    print "FWHM_Si=", sp.specfit.measure_approximate_fwhm()
    C_fit=line_fitter("C", i,guesses_C, limits_C, limited_C, tied_C, xmin_C, xmax_C,magorder=-1*mag_order)
    print "FWHM_C=", sp.specfit.measure_approximate_fwhm()
    xmin_fe=4420
    xmax_fe=4700
    template_file="./templates/fe2_Op.dat"
    fe_template = pyspeckit.Spectrum(template_file)
    fe_template.xarr.units='angstroms'
    fe_template.xarr.xtype='wavelength'
    lambda_fe=4700
    index_sp=np.argmin( np.abs( lambda_fe - sp.xarr) )
    index_fe=np.argmin( np.abs( lambda_fe - fe_template.xarr) )
    flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
    fe_template.data*=flux/fe_template.data[index_fe]
    fe_Hb=blue_bump_fitter("optical",i,xmin_fe,xmax_fe,fe_template, balmer_template)
    
    #print limits_Hbeta
    Hbeta_fit=line_fitter("Hbeta", i, guesses_Hbeta, limits_Hbeta, limited_Hbeta, tied_Hbeta, xmin_Hbeta, xmax_Hbeta,offset=-0.05, magorder=-1*mag_order)
    print "FWHM_Hbeta=", sp.specfit.measure_approximate_fwhm()
    Halpha_fit=line_fitter("Halpha", i, guesses_Halpha, limits_Halpha, limited_Halpha, tied_Halpha, xmin_Halpha, xmax_Halpha,offset=-0.05, magorder=-1*mag_order)
    print "FWHM_Halpha=", sp.specfit.measure_approximate_fwhm()
    print len(Halpha_fit)
    print len(Hbeta_fit)
    print len(C_fit)
    print len(Mg_fit)
    print len(fe_fit)
    print len(balmer_template.data)
    total=continuous +  Halpha_fit + Hbeta_fit + C_fit + fe_fit  + balmer_data + Si_fit + Mg_fit
    """
    fitmin=sp.specfit.xmin
    fitmax=sp.specfit.xmax
    pylab.figure()
    pylab.xlim(xmin_MgII,xmax_MgII)
    pylab.plot(backup.xarr,backup.data)
    pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals + fe_fit[fitmin:fitmax] + mg_fit[fitmin:fitmax] )
    pylab.plot(sp.xarr,sp.data - mg_fit )
    
    pylab.savefig("Hbeta_fe_fit_"+str(i)+".pdf")
    """
    pylab.figure()
    pylab.ylim(ymin=-1.3*np.abs(copy1.data.min()),ymax=1.1*copy1.data.max())
    pylab.xlim(xmin=1230,xmax=9900)
    pylab.ylabel(r'$10^{-'+str(mag_order)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.xscale('log')
    pylab.plot(sp.xarr,copy3.data,'k')
    pylab.plot(sp.xarr, total,'r')
    total=continuous + fe_fit + balmer_data 
    pylab.plot(sp.xarr, total,'b')
    total=continuous + balmer_data 
    pylab.plot(sp.xarr, total,'g')
    total=continuous
    pylab.plot(sp.xarr, total,'gray')
    
    plot_file="./plots/superposition_" + str(i) + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)






    print L_cont



"""

sp.measure(z = 1.55, fluxnorm = 1e-16)


y = sp.plotter.ymax * 0.85    # Location of annotations in y

for i, line in enumerate(sp.measurements.lines.keys()):

# If this line is not in our database of lines, don't try to annotate it
    if line not in sp.speclines.optical.lines.keys(): continue

    x = sp.measurements.lines[line]['modelpars'][1]   # Location of the emission line


# Print out spectral line information
print "Line   Flux (erg/s/cm^2)     Amplitude (erg/s/cm^2)    FWHM (Angstrom)   Luminosity (erg/s)"
for line in sp.measurements.lines.keys():
    print line, sp.measurements.lines[line]['flux'], sp.measurements.lines[line]['amp'], sp.measurements.lines[line]['fwhm'], \
        sp.measurements.lines[line]['lum']

# Had we not supplied the objects redshift (or distance), the line
# luminosities would not have been measured, but integrated fluxes would
# still be derived.  Also, the measurements class separates the broad and
# narrow H-alpha components, and identifies which lines are which. How nice!

sp.specfit.plot_fit()
"""



