from matplotlib import pylab
import numpy as np
import pyspeckit
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from scipy.ndimage.filters import gaussian_filter as gauss_conv
from astropy.cosmology import FlatLambdaCDM
from fitcode.line import *
execfile("./constraints.cfg") 

#from pyspeckit.spectrum.models.template import template_fitter as temp_fitter
spec_division=["OP","UV","FUV"] # Do not change the order, The Far UV (FUV) fit 
                                # depens on the UV Fe template fit. 
xmin=1250
xmax=9000
n_spec=1
sp_to_run=[1,3,4,8,9,11,15,16,18,20,22,24,27]
#sp_to_run=[15,16,18,20,22,24,27]
##sp_to_run=[1,3,4,8,9,11,12,15,16,18,20,22,24,27]
sp_to_run=[1]

pylab.ion()



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
            continuous,wlmin_UV=continuous_substraction( i, sp, mag_order,UV_limits,w)
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
            continuous,wlmin_FUV=continuous_substraction( i, sp, mag_order,FUV_limits,w)
            
        if w=="OP": continuous,wlmin_OP=continuous_substraction( i, sp, mag_order,OP_limits,w)
        #-------------continuous subtraction-------------------#
                

        #------------------- -------Shared features for fitting. Finished ---------------------------#

        

        
        
        # -------------------------------------------------Halpha plus Hbeta--------------------------------- 
        if w=="OP": 
            continuous_OP=continuous
            sp.data=sp.data - continuous_OP 
            xmin_fe=4600
            xmax_fe=4800
            template_file=fe_template_OP
            fe_template = pyspeckit.Spectrum(template_file)
            xmin_fe_template=fe_template.xarr[0]
            fe_template.xarr.units='angstroms'
            fe_template.xarr.xtype='wavelength'
            lambda_fe=4500
            index_sp=np.argmin( np.abs( lambda_fe - sp.xarr) )
            index_fe=np.argmin( np.abs( lambda_fe - fe_template.xarr) )
            flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
            fe_template.data*=np.abs(flux/fe_template.data[index_fe])
            
            pylab.rcParams["figure.figsize"]=16,6
            
            xmin_fe=fe_template.xarr[0]
            xmax_fe=fe_template.xarr[-1]
            print "xmin,xmax=   ",xmin_fe,xmax_fe
        
            fe_Hb=blue_bump_fitter(sp,"optical",i,xmin_fe,xmax_fe,fe_template,mag_order,do_fit=fit_fe_OP)
            pylab.rcParams["figure.figsize"]=8,6
            
            Hbeta_fit=line_fitter(sp, "Hbeta", i, guesses_Hbeta, limits_Hbeta, limited_Hbeta, tied_Hbeta, xmin_Hbeta, xmax_Hbeta,offset=-0.05, magorder=mag_order, do_fit=fit_Hbeta)
            if fit_Hbeta: print "FWHM_Hbeta=", sp.specfit.measure_approximate_fwhm()
            Halpha_fit=line_fitter(sp,"Halpha", i, guesses_Halpha, limits_Halpha, limited_Halpha, tied_Halpha, xmin_Halpha, xmax_Halpha,offset=-0.05, magorder=mag_order, do_fit=fit_Halpha)
            if fit_Halpha: print "FWHM_Halpha=", sp.specfit.measure_approximate_fwhm()
            
        #-------------------------------------------------Halpha plus Hbeta--------------------------------- 
             
    
    

        # -------------------------------------------------MgII plus FeII UV plus Balmer Cont--------------------------------- 
        if w=="UV":
            continuous_UV=continuous
            sp.data=sp.data - continuous_UV 
            #---------Normalizing Fe Template -------------
            xmin_bb=3640
            xmax_bb=3650
            lambda_0=3675
            index=np.argmin( np.abs( lambda_0 - sp.xarr) )
            flux=np.mean(sp.data[index -30: index + 30])
            template_file=fe_template_UV
            fe_template_UV = pyspeckit.Spectrum(template_file)
            fe_template_UV.xarr.units='angstroms'
            fe_template_UV.xarr.xtype='wavelength'
            #fe_template_UV.crop(1555.0,fe_template_UV.xarr[-1],units="angstroms")
            lambda_fe=2700
            index_sp=np.argmin( np.abs( lambda_fe - sp.xarr) )
            index_fe=np.argmin( np.abs( lambda_fe - fe_template_UV.xarr) )
            flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
            fe_template_UV.data*=flux/fe_template_UV.data[index_fe]
            #---------Normalizing Fe Template -------------

            
        
            #-------------Fitting the Blue Bump------------------------
            xmin_fe=2200
            xmax_fe=3500
            
            fe_fit=blue_bump_fitter(sp,"SBB",i,xmin_fe,xmax_fe,fe_template_UV,mag_order, do_fit=fit_fe_UV)
            pylab.rcParams["figure.figsize"]=8,6            
            Mg_fit=line_fitter(sp, "MgII", i, guesses_MgII, limits_MgII, limited_MgII, tied_MgII, xmin_MgII, xmax_MgII,offset=-0.05, magorder=mag_order,do_fit=fit_Mg)
            if fit_Mg: print "FWHM_Mg=", sp.specfit.measure_approximate_fwhm()
            FWHM_Mg=sp.specfit.measure_approximate_fwhm()
            FWHMV_Mg=wl_to_kms(FWHM_Mg,2800)
            #-------------Fitting the Blue Bump------------------------


            
            #-------------Smoothing Balmer template with MgII width ------------------------
            smooth_f=FWHMV_Mg/600.0 
            sigmav_Mg=FWHMV_Mg/(2.355)
            sigma=kms_to_wl(np.sqrt(sigmav_Mg*sigmav_Mg-600.0*600.0),3800) 
            balmer_tot=gauss_conv(balmer_tot, sigma, order=0, output=None, mode='reflect', cval=0.0)
            print "smooth", smooth_f,sigmav_Mg
            balmer_data=balmer_tot
            #-------------Smoothing Balmer template with MgII width ------------------------

            
            #line_fitter(sp, "HeI", i, guesses_HeI, limits_HeI, limited_HeI, tied_HeI, xmin_HeI, xmax_HeI)
            
            
            
            
        if w=="FUV":
            continuous_FUV=continuous_UV
            sp.data=sp.data - continuous_FUV 
            sp.data=sp.data - fe_fit
            Si_fit=line_fitter(sp,"SiIV", i,guesses_SiIV, limits_SiIV, limited_SiIV, tied_SiIV, xmin_SiIV, xmax_SiIV, magorder=mag_order,do_fit=fit_Si)
            if fit_Si: print "FWHM_Si=", sp.specfit.measure_approximate_fwhm()
            C_fit=Si_fit*0.0
            C_fit=line_fitter(sp, "C", i,guesses_C, limits_C, limited_C, tied_C, xmin_C, xmax_C,magorder=mag_order, do_fit=fit_C)
            if fit_C: print "FWHM_C=", sp.specfit.measure_approximate_fwhm()
            
    number,name,mag,group,redshift=np.genfromtxt("quasar_data.txt", dtype="|S10",unpack=True)       
    del(number,mag,redshift)
            
    sp = pyspeckit.Spectrum(spectrum_file)
    
    # -----------set up units properly------------------#        
    mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    sp.xarr.units='angstroms'
    sp.xarr.xtype = 'wavelength'
    sp.units = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $cm^{-1}$'
    sp.data /= 10**(mag_order)    
    #-------------- set up units properly------------#
    
    xmin_fe_template=fe_template=fe_template.xarr[0]
    
    copy1=sp.copy()
    continuous=0.0*sp.data
    arg_UV=np.argmin( np.abs(sp.xarr-wlmin_UV) )
    arg_OP=np.argmin( np.abs(sp.xarr - wlmin_OP) )
    arg_OP=np.argmin( np.abs(sp.xarr - xmin_fe_template)  )
    #continuous[:arg_UV]=continuous_FUV[:arg_UV]
    #continuous[arg_UV:arg_OP]=continuous_UV[arg_UV:arg_OP]
    continuous[:arg_OP]=continuous_UV[:arg_OP]
    continuous[arg_OP:]=continuous_OP[arg_OP:]
    
    
    #----Superposition plot --------#
    total=continuous +  Halpha_fit + Hbeta_fit + C_fit + fe_fit  + fe_Hb + balmer_data + Si_fit + Mg_fit
    
    pylab.rcParams["figure.figsize"]=16,6
        
    copy1.crop(1400,1600)
    pylab.figure()
    pylab.ylim(ymin=-1.0*np.abs(copy1.data.min()),ymax=1.1*copy1.data.max())
    pylab.xlim(xmin=wlmin,xmax=wlmax)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.yscale('linear')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(sp.xarr, total,'r')
    total=continuous + fe_fit  + balmer_data + fe_Hb
    pylab.plot(sp.xarr, total,'b')
    total=continuous + balmer_data 
    pylab.plot(sp.xarr, total,'g')
    total=continuous
    pylab.plot(sp.xarr, total,'gray')
    
    plot_file="./plots/superposition_" + group[i]+"_"+ name[i] + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)
    pylab.show()
    #----Superposition plot --------#



    
    #---- FUV zoom --------#
    
    pylab.figure()
    xmin=sp.xarr.min()
    xmax=2300
    copy1=sp.copy()
    copy1.crop(sp.xarr.min(),2300)

    pylab.ylim(ymin=0,ymax=1.05*copy1.data.max())

    pylab.xlim(xmin=xmin,xmax=xmax)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.yscale('linear')
    pylab.plot(sp.xarr,sp.data,'k')
    total=continuous  + C_fit + fe_fit + balmer_data + Si_fit     
    pylab.plot(sp.xarr, total,'r')
    total=continuous + balmer_data +fe_fit  
    pylab.plot(sp.xarr, total,'b')    

    
    total=continuous + balmer_data 
    pylab.plot(sp.xarr, total,'g')
    
    total=continuous
    pylab.plot(sp.xarr, total,'gray')
    
    plot_file="./plots/FUV_zoom_" + group[i]+"_"+ name[i] + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)

    pylab.rcParams["figure.figsize"]=16,6
    pylab.show()
    
    #----FUV zoom --------#



    #----UV zoom --------#
    pylab.figure()
    xmin=2200
    xmax=4050
    copy1=sp.copy()
    copy1.crop(2600,4000)

    pylab.ylim(ymin=0,ymax=1.05*copy1.data.max())

    pylab.xlim(xmin=xmin,xmax=xmax)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.yscale('linear')
    pylab.plot(sp.xarr,sp.data,'k')
    pylab.plot(sp.xarr, total,'r')
    total=continuous + fe_fit  + balmer_tot
    pylab.plot(sp.xarr, total,'b')
    
    
    total=continuous + balmer_data 
    pylab.plot(sp.xarr, total,'g')
    
    total=continuous
    pylab.plot(sp.xarr, total,'gray')
    
    plot_file="./plots/UV_zoom_" + group[i]+"_"+ name[i] + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)

    pylab.rcParams["figure.figsize"]=16,6
    pylab.show()
    
    #----UV zoom --------#
    

    #----Optical zoom --------#
    
    pylab.figure()
    xmin=4000
    xmax=7000
    copy1=sp.copy()
    copy1.crop(4000,7000)

    pylab.ylim(ymin=0,ymax=1.05*copy1.data.max())

    pylab.xlim(xmin=xmin,xmax=xmax)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.yscale('linear')
    pylab.plot(sp.xarr,sp.data,'k')
    total=continuous +  Halpha_fit + Hbeta_fit +  fe_Hb 
    pylab.plot(sp.xarr, total,'r')
    total=continuous + fe_Hb  
    pylab.plot(sp.xarr, total,'b')    
    total=continuous
    pylab.plot(sp.xarr, total,'gray')
    
    plot_file="./plots/OP_zoom_" + group[i]+"_"+ name[i] + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)

    pylab.rcParams["figure.figsize"]=16,6
    pylab.show()
    
    #----Optical zoom --------#
    
    pylab.figure()
    pylab.rcParams["figure.figsize"]=16,6
    #----residual plot --------#
    total=continuous +  Halpha_fit + Hbeta_fit + C_fit + fe_fit  + fe_Hb + balmer_data + Si_fit + Mg_fit
    copy1=sp.copy()
    copy1.crop(1400,1600)
    pylab.ylim(ymin=-1.3*copy1.data.min(),ymax=1.0*copy1.data.max())
    pylab.xlim(xmin=sp.xarr[0],xmax=sp.xarr[-1])
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.xscale('log')
    pylab.plot(sp.xarr ,sp.data - total,'k')
    plot_file="./plots/no_continuous_" + group[i]+"_"+ name[i] + ".png"
    
    pylab.savefig(plot_file)
    

    


