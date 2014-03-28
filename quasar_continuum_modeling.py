from matplotlib import pylab
import numpy as np
import pyspeckit
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from scipy.ndimage.filters import gaussian_filter as gauss_conv
from astropy.cosmology import FlatLambdaCDM
cosmology=FlatLambdaCDM(H0=70,Om0=0.3)
from astropy import units as u
from fitcode.line import *

#------Executing the configuration file that includes the requiered parameters---#
execfile("./constraints.cfg") 
#------Executing the configuration file that includes the requiered parameters---#


#from pyspeckit.spectrum.models.template import template_fitter as temp_fitter

Mpc=3.08567758e24 #centimeters

mAB_nuv_galex=np.array(["19.789","19.171","nan  ","18.72","17.657","19.452","18.914","20.735",
                    "nan   ","19.337","21.511","20.367","19.585","19.01", "nan   ","20.261",
                    "21.787","20.743","20.809","21.394","20.473","20.405","nan   ","21.502",
                    "20.603","21.059","nan   ","22.837","21.049","21.292"])


mAB_fuv_galex=np.array(["nan", "nan", "nan", "nan","18.689","20.633","nan","20.759","nan",
                        "20.53","21.53","22.473","21.422","nan","nan","nan","23.171","22.04",
                        "nan","nan","nan","nan","nan","nan","nan","nan","nan","nan","20.57",
                        "nan"])









20.57


color_index=np.array([0.02705,0.02097,0.04293,0.02034,0.01067,0.04108,0.02407,0.03662,0.05813,0.08435,
             0.03108,0.0356,0.02471,0.02175,0.08482,0.03347,0.04569,0.03438,0.0332,0.04171,
             0.02575,0.02631,0.03489,0.0241,0.04068,0.0279,0.0867,0.13133,0.04065,0.12069])

R_nuv=3.0*3.1#7.347#6.738
R_fuv=6.892

A_nuv=R_nuv*color_index

xmin=1250
xmax=9000

mAB_nuv_galex=np.array([float(x) for x in mAB_nuv_galex])
mAB_fuv_galex=np.array([float(x) for x in mAB_fuv_galex])
mAB_nuv_galex0=mAB_nuv_galex
mAB_nuv_galex=mAB_nuv_galex - A_nuv
mAB_fuv_galex=mAB_fuv_galex - A_nuv


#pylab.ion()


for i in sp_to_run:

    spectrum_file="./AGN_spectra/spectrum" + str(i) + ".txt"
    #-------------continuous fitting-------------------#
    exclude_file = "./excludes/exclude_cont.txt"
    exclude_cont=np.loadtxt(exclude_file,skiprows=2)
    sp = pyspeckit.Spectrum(spectrum_file)
    
            
    k=1
    
    for w in spec_division:
    


        print("spectrum number , region = " + str(i) + " ,  " + w)
        sp = pyspeckit.Spectrum(spectrum_file)
    
        #------------opening spectrum i---------------------#
    
        
        # -----------set up units properly------------------#        
        mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
        sp.xarr.units='angstroms'
        sp.xarr.xtype = 'wavelength'
        sp.units = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $cm^{-1}$'
        sp.data /= 10**(mag_order)    
        #-------------- set up units properly------------#
        
        
        wlmin=sp.xarr[0]
        wlmax=sp.xarr[-1]
    
    
        pylab.yscale('log')
        pylab.xscale('log')
        
        #pylab.xscale('log',subsx=[3,4])
        pylab.rcParams["figure.figsize"]=16,6
        sp.plotter(figure=1,xmin=wlmin,xmax=wlmax,ymin=1.1*sp.data.min(),ymax=3.5*sp.data.max())
        pylab.close()
        
        #-------------continuous subtraction-------------------#
        if w=="UV": 
            continuous,wlmin_UV,model_contUV=continuous_substraction( i, sp, mag_order,UV_limits,w)
            backup=sp.copy()
            sp.data=sp.data - continuous
            balmer_template,balmer_tot, index=balmer_normalization(sp,balmer_cont_template,balmer_lines_template)
            
            sp.data=backup.data
            
            cont=sp.copy()
            cont.data=continuous
            
            cont,balmer_template,scale_balmer=total_continuous_fit(i,sp,cont,balmer_template,mag_order,w)
            continuous=cont.data
            balmer_tot=scale_balmer*balmer_tot
            k=0
            
            #-------extracting balmer continuum (not lines yet)-------#
            sp.data[:index]=sp.data[:index] - balmer_template.data[:index] 
            #-------extracting balmer continuum (not lines yet)-------#

        if w=="FUV": 
            sp.data[:index]=sp.data[:index] - balmer_template.data[:index] 
            continuous,wlmin_FUV,model_contFUV=continuous_substraction( i, sp, mag_order,FUV_limits,w)
            
        if w=="OP": continuous,wlmin_OP,model_contOP=continuous_substraction( i, sp, mag_order,OP_limits,w)
        #-------------continuous subtraction-------------------#
                

        #------------------- -------Shared features for fitting. Finished ---------------------------
        
        

        if w=="OP": 
            continuous_OP=continuous
            sp.data=sp.data - continuous_OP 
            
    


        if w=="UV":
            continuous_UV=continuous
            
            sp.data=sp.data - continuous_UV 
            
            
            
        if w=="FUV":
            continuous_FUV=continuous_UV
            sp.data=sp.data - continuous_FUV 
            

    number,name,mag,group,redshift=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
    
    
    
    
    redshift=np.array([float(x) for x in redshift])
    mag=np.array([float(x) for x in mag])
    dl=a=cosmology.luminosity_distance(redshift).to(u.cm)
    
    F_nuv_galex=galex_mAB_Fnuv(mAB_nuv_galex)*u.erg/(u.second*u.cm*u.cm*u.angstrom)
    galex_central_wl=2315.7*u.angstrom/(1+redshift)
    L_wl=((1+redshift)*4*np.pi*F_nuv_galex*dl.to(u.cm)*dl.to(u.cm)).to(u.erg/(u.cm*u.s))
    L_nu=(L_wl*np.power(2315.7*u.angstrom/(1+redshift),2)/(3e10*u.cm/u.second)).to(u.erg)
    




    L_wl=np.multiply( L_wl ,10**(-1*mag_order) )
    
    L_wl=L_wl.value
    galex_central_wl=galex_central_wl.value



    F_fuv_galex=galex_mAB_Ffuv(mAB_fuv_galex)*u.erg/(u.second*u.cm*u.cm*u.angstrom)
    galex_central_wlf=1538.6*u.angstrom/(1+redshift)
    L_wlf=((1+redshift)*4*np.pi*F_fuv_galex*dl.to(u.cm)*dl.to(u.cm)).to(u.erg/(u.cm*u.s))
    L_nuf=(L_wlf*np.power(2315.7*u.angstrom/(1+redshift),2)/(3e10*u.cm/u.second)).to(u.erg)
    L_wlf=np.multiply( L_wlf ,10**(-1*mag_order) )    
    L_wlf=L_wlf.value
    galex_central_wlf=galex_central_wlf.value
    




    print L_wl[i], galex_central_wl[i]
    

    model=np.loadtxt("./model/cont_model_" + str(i) + ".txt")
                                                                       
    #Lmodel=fdata=np.interp(sp.xarr,model[:,0],model[:,1])             
                                                                       
    #Fmodel=np.multiply( Lmodel ,10**(-1*magorder) )                   
    model[:,1]=np.multiply( model[:,1] ,10**(-1*mag_order) )    
    
    sp = pyspeckit.Spectrum(spectrum_file)
    
    # -----------set up units properly------------------#        
    mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    sp.xarr.units='angstroms'
    sp.xarr.xtype = 'wavelength'
    sp.units = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $cm^{-1}$'
    sp.data /= 10**(mag_order)    
    #-------------- set up units properly------------#
    
    print "hola 1"
    
    copy1=sp.copy()
    continuous=0.0*sp.data
    arg_UV=np.argmin( np.abs(sp.xarr-wlmin_UV) )
    #arg_OP=np.argmin( np.abs(sp.xarr - wlmin_OP) )
    arg_OP=np.argmin( np.abs(sp.xarr - 4050.0)  )
    #continuous[:arg_UV]=continuous_FUV[:arg_UV]
    #continuous[arg_UV:arg_OP]=continuous_UV[arg_UV:arg_OP]
    continuous[:arg_OP]=continuous_UV[:arg_OP]
    continuous[arg_OP:]=continuous_OP[arg_OP:]
    
    print "hola 2"
    
    #----Superposition plot --------#
    total=(continuous +  balmer_tot)
    
    pylab.rcParams["figure.figsize"]=16,6
        
    copy1.crop(1400,1600)
    pylab.figure()
    #pylab.ylim(ymin=-1.0*np.abs(copy1.data.min()),ymax=1.1*copy1.data.max())
    #pylab.xlim(xmin=wlmin,xmax=wlmax)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.yscale('linear')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(model_contFUV[:,0], model_contFUV[:,1],'r')
    print "hola 3"
    pylab.plot(galex_central_wl[i], L_wl[i],'bo',markersize=5.0)
    pylab.plot(galex_central_wlf[i], L_wlf[i],'bo',markersize=5.0)
    print "hola 4"
    #pylab.plot(model[:,0], model[:,1],'g')
    plot_file="./plots/GALEX_" + group[i]+"_"+ name[i] + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)
    print "hola 5"
    pylab.show()
    print "hola 6"
    #----Superposition plot --------#



    


