from matplotlib import pylab
import numpy as np 
import pyspeckit
import os
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from scipy.ndimage.filters import gaussian_filter as gauss_conv
#import quasar_fit_continuum

def continuous_substraction( spec_number, sp, magorder,cont_limits):
    
    model=np.loadtxt("./model/cont_model_" + str(spec_number) + ".txt")

    Lmodel=fdata=np.interp(sp.xarr,model[:,0],model[:,1])

    
    

    Fmodel=np.multiply( Lmodel ,10**(-1*magorder) )
    #printA Lmodel,Fmodel, sp.data, magorder

    exclude_file = "./excludes/exclude_cont.txt"
    exclude_cont=np.loadtxt(exclude_file,skiprows=2)
    
    cont=exclude_cont[spec_number,:][1:]
    print "cont=  ", cont
    arg_min=np.argmin( np.abs( cont - cont_limits[0] ) )
    if cont[arg_min] - cont_limits[0] > 0 : arg_min= arg_min + 1
    arg_max=np.argmin( np.abs( cont - cont_limits[1] ) )
    if cont[arg_max] - cont_limits[0] < 0 : arg_max= arg_max - 1
    cont=cont[arg_min:arg_max+1]
    
    print "cont=  ", cont
    
    Nregions=np.int( len(cont)/2 )
    
    print "Nreg= ", Nregions

    arg_list=np.ones(len(cont),dtype="int")
    i=0

    for wl in cont:
        arg= np.argmin(np.abs(sp.xarr - wl) ) 
        arg_list[i]=arg
        i=i+1
    ratio=np.zeros(Nregions)
    
    for i in range(Nregions):
        ratio[i]=np.abs( np.mean(sp.data[ arg_list[2*i]:arg_list[2*i + 1] ])/np.mean(Fmodel[ arg_list[2*i]:arg_list[2*i + 1] ]) )
        
        mean_ratio=np.mean(ratio)
    
    
    Fmodel*=mean_ratio
    
    pylab.figure()
    pylab.ylabel(r'$10^{'+str(magorder)+'}$ erg s$^{-1}$ $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(sp.xarr,Fmodel,'r')
    pylab.show()

    sp.data=sp.data - Fmodel 
    return Fmodel, cont[arg_min]

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

def wl_to_kms(wl,line_center):
    c=3e5
    kms=wl*c/line_center
    return kms


    
def blue_bump_fitter(sp,line_name,spec_number,xmin,xmax,fe_template, magorder):
    
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
            return new_fe_template


        


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
        
        sp.specfit.fullsizemodel()  
        
        fe_template.data*=sp.specfit.parinfo.values[0]
        fe_template.xarr+=sp.specfit.parinfo.values[1]
        new_template=fe_template.copy()
        
        def suavizado(xarr, smooth_factor_fe,templat=new_template.copy()):#,scale_balmer):,fe_template=fe_template.copy(), balmer_template=balmer_template.copy()):
            
            
            data=smooth(templat.data,smooth_factor_fe , downsample=False)   
            templat.data=data
            n_template =pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                        templat.xarr.as_unit('angstroms'), 
                                                        templat.data)
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
    pylab.ylabel(r'$10^{'+str(magorder)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
    pylab.savefig(plot_file)
    # -----------------  ploting fit---------------------
    
    
    # ----------plotting residuals-------------
    
    pylab.figure()
    pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals[fitmin:fitmax].max())
    pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
    pylab.ylabel(r'$10^{-'+str(magorder)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    
    pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
    plot_file="./plots/" + line_name +"_res_" + str(spec_number) + ".png"
    
    pylab.savefig(plot_file)
    # ----------plotting residuals-------------#
    
    #----------redifining spectrum ------------------#
    sp.data=sp.data - sp.specfit.fullmodel
    #----------redifining spectrum ------------------#
    return sp.specfit.fullmodel
        
def line_fitter(sp,line_name, spec_number,guesses, limits, limited, tied, xmin, xmax, offset=-0.4, magorder=-16):
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
    sp.measure(z = 0.00, fluxnorm = 10**(magorder))
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
    pylab.ylabel(r'$10^{'+str(magorder)+'}$ erg s$^{-1}$  $cm^{-1}$')
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
        
        print_line = '{0}  {1}  {2}  {3} {4} {5}\n'.format(line, sp.measurements.lines[line]['flux'],sp.measurements.lines[line]['amp'],sp.measurements.lines[line]['fwhm'],sp.measurements.lines[line]['lum'], sp.measurements.lines[line]['modelpars'][1] )
        f.write(print_line)
    f.close()
    return sp.specfit.fullmodel
    #---------------measuring----------------------#


def fe_fitter(sp,line_name,spec_number,xmin,xmax,template,magorder):
    
    #-----------------------------defining fe template fitter ------------------------------------#
    #template_file="./templates/fe2_UV.dat"
    #template = pyspeckit.Spectrum(template_file)

    for i in range(5):
        
        def fe_model(xarr, scale_fe,shift_fe, smooth_factor_fe):

            

            new_fe_template = scale_fe*pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                                       template.xarr.as_unit('angstroms')+shift_fe, 
                                                                       template.data )
        
            
            

            new_template=  new_fe_template 
            return new_fe_template


        


        feII = pyspeckit.spectrum.models.model.SpectralModel(fe_model,
                                                                   3, 
                                                                   parnames=['scale_fe','shift_fe','smooth_fe'],#,'scale_balmer'], 
                                                                   parlimited=[(True,True),(False,False),(True,False)],#,(False,False)],  
                                                                   parlimits=[(0.1,10.0), (-10.0,10.0),(0,0)],#,(0.9,1.1)], 
                                                                   shortvarnames=(r'Af',r'\Delta xf',r'sigma_xf'),#,r'Ab'),
                                                                   centroid_par='shift_fe' )
        feII.__name__ = "fe_template"



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
            
            
        
        sp.Registry.add_fitter('fe_template',feII,3,multisingle='multi')
        sp.specfit(fittype='fe_template',exclude=exclude,guesses=[1.0,0,1.0],xmin=xmin,xmax=xmax)
        
        sp.specfit.fullsizemodel()  
        
        template.data*=sp.specfit.parinfo.values[0]
        template.xarr+=sp.specfit.parinfo.values[1]
        new_template=template.copy()
        
        def suavizado(xarr, smooth_factor_fe,template=new_template.copy()):
            
            
            data=smooth(template.data,smooth_factor_fe , downsample=False)   
            template.data=data
            n_template =pyspeckit.interpolation._interp(xarr.as_unit('angstroms'), 
                                                        template.xarr.as_unit('angstroms'), 
                                                        template.data)
            template.data=data
            
            
            
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
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(sp.xarr,sp.specfit.fullmodel,'r')
    pylab.savefig(plot_file)
    # -----------------  ploting fit---------------------
    
    
    # ----------plotting residuals-------------
    
    pylab.figure()
    pylab.ylim(ymin=-1.2*np.abs(sp.data[fitmin:fitmax].min()),ymax=1.1*sp.specfit.residuals[fitmin:fitmax].max())
    pylab.xlim(xmin=sp.xarr[fitmin],xmax=sp.xarr[fitmax])
    pylab.ylabel(r'$10^{-'+str(magorder)+'}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    pylab.xlabel(r'$\AA$')
    
    pylab.plot(sp.xarr[fitmin:fitmax],sp.specfit.residuals[fitmin:fitmax],'k')
    plot_file="./plots/" + line_name +"_res_" + str(spec_number) + ".png"
    
    pylab.savefig(plot_file)
    # ----------plotting residuals-------------#

    wlmin=np.amin(np.abs(sp.xarr-xmin))
    wlmax=np.amin(np.abs(sp.xarr-xmin))
    model=sp.specfit.fullmodel
    model[:xmin+1]=0.0
    model[xmax:]=0.0
    
    #----------redifining spectrum ------------------#
    sp.data=sp.data - model
    #----------redifining spectrum ------------------#
    return model

