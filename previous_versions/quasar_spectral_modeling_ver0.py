
from matplotlib import pylab
#import numpy as np
import pyspeckit
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpolation
from scipy.ndimage.filters import gaussian_filter as gauss_conv
from astropy.cosmology import FlatLambdaCDM
from fitcode.line import *
from fitcode.constraints import *

#from pyspeckit.spectrum.models.template import template_fitter as temp_fitter
xmin=1250
xmax=9000
n_spec=1
#sp_to_run=[1]
#sp_to_run=[1,3,4,8,9,11,12,15,16,18,20,22,24,27]
sp_to_run=[3,4,8,9,11,15,16,18,20,22,24,27]
#sp_to_run=[15,16,18,20,22,24,27]
#sp_to_run=[15,16,18,20,22,24,27]
#sp_to_run=[1,17,20]
#sp_to_run=[1]

pylab.ion()



for i in sp_to_run:
    
    pylab.rcParams["figure.figsize"]=16,6
    #------------------- -------Shared features for fitting. Finished -----------------------#
    #------------opening spectrum i---------------------#
    spectrum_file='../AGN_spectra/spectrum' + str(i) + '.txt'
    sp = pyspeckit.Spectrum(spectrum_file)
    #sp.smooth(2)
    #------------opening spectrum i---------------------#
    
    
    # -----------set up units properly------------------#
    mag_order=np.int((1)*np.round(np.log10(np.mean(sp.data))))
    sp.xarr.units='angstroms'
    sp.xarr.xtype = 'wavelength'
    sp.units = r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$ $cm^{-1}$'
    sp.data /= 10**(mag_order)    
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

    #limit_file = "./cont_limits.txt"
    #cont_limits=np.loadtxt(limit_file)
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

    continuous,L_cont=continuous_substraction( i, sp, mag_order,UV_limits)
    copy1.crop(1400,1600)
    pylab.ylim(ymin=-1.3*np.abs(copy1.data.min()),ymax=1.1*copy1.data.max())
    pylab.xlim(xmin=xmin,xmax=6900)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.xscale('log')
    pylab.plot(sp.xarr,sp.data,'k')
    plot_file="./plots/no_continuous_" + str(i) + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)
    
    #-------------continuous fitting-------------------#

    
    #------------------- -------Shared features for fitting. Finished ---------------------------#




    fe_Hb=0.0*sp.data
    """ -------------------------------------------------Halpha plus Hbeta--------------------------------- 
    xmin_fe=4600
    xmax_fe=4800
    template_file="./templates/fe2_Op.dat"
    fe_template = pyspeckit.Spectrum(template_file)
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
    
    fe_Hb=fe_fitter("optical",i,xmin_fe,xmax_fe,fe_template,mag_order)
    pylab.rcParams["figure.figsize"]=8,6
    #print limits_Hbeta
    Hbeta_fit=line_fitter(sp, "Hbeta", i, guesses_Hbeta, limits_Hbeta, limited_Hbeta, tied_Hbeta, xmin_Hbeta, xmax_Hbeta,offset=-0.05, magorder=mag_order)
    print "FWHM_Hbeta=", sp.specfit.measure_approximate_fwhm()
    Halpha_fit=line_fitter(sp,"Halpha", i, guesses_Halpha, limits_Halpha, limited_Halpha, tied_Halpha, xmin_Halpha, xmax_Halpha,offset=-0.05, magorder=mag_order)
    print "FWHM_Halpha=", sp.specfit.measure_approximate_fwhm()

    -------------------------------------------------Halpha plus Hbeta--------------------------------- 
    """ 
    
    


    
    xmin_CIII=1940
    xmax_CIII=2010

    #fe_fitter(sp,"C",i,xmin_CIII,xmax_CIII)
    #line_fitter(sp, "C", i,guesses_C, limits_C, limited_C, tied_C, xmin_C, xmax_C,magorde=-1*mag_order)
    

    #line_fitter(sp, "HeI", i, guesses_HeI, limits_HeI, limited_HeI, tied_HeI, xmin_HeI, xmax_HeI)
    
    
    xmin_bb=3640
    xmax_bb=3650
    lambda_0=3675
    index=np.argmin( np.abs( lambda_0 - sp.xarr) )
    flux=np.mean(sp.data[index -30: index + 30])
    
    
    #sp.data=sp.data - sp.specfit.fullmodel
    template_file="./templates/fe2_UV2.dat"
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
    balmer_lines_file="./templates/balmer_lines2.txt"
    balmer_template = np.loadtxt(template_file)
    balmer_lines=np.loadtxt(balmer_lines_file)
    balmer_data=np.interp(sp.xarr,balmer_template[:,0], balmer_template[:,1])
    balmer_lines_data=np.interp(sp.xarr,balmer_lines[:,0], balmer_lines[:,1])
    balmer_template=pyspeckit.Spectrum(data=balmer_data,xarr=sp.xarr)
    balmer_template.xarr.units='angstroms'
    balmer_template.xarr.xtype='wavelength'


    balmer_lines=pyspeckit.Spectrum(data=balmer_lines_data,xarr=sp.xarr)
    balmer_lines.xarr.units='angstroms'
    balmer_lines.xarr.xtype='wavelength'

    
    lambda_balmer=3645
    index=np.argmin( np.abs( lambda_balmer - sp.xarr) )
    index_hi=np.argmin( np.abs( 4000.0 - sp.xarr) )
    flux=np.mean(sp.data[ index - 30 : index + 30 ])
    balmer_lines.data*=flux/balmer_template.data[-1]
    balmer_template.data*=flux/balmer_template.data[-1]

    balmer_tot=np.zeros(len(balmer_template.data))
    balmer_tot[:index]=balmer_template.data[:index]
    balmer_tot[index:index_hi]=balmer_lines.data[index:index_hi]
    balmer_tot1=balmer_tot
    
    backup=sp.copy()
    sp.data[:index]=sp.data[:index] - balmer_template.data[:index]
    
    xmin_fe=2200
    xmax_fe=3500

    backup=sp.copy()
    #fe_fit=0.0*Hbeta_fit
    fe_fit=blue_bump_fitter(sp,"SBB",i,xmin_fe,xmax_fe,fe_template,mag_order)

    pylab.rcParams["figure.figsize"]=8,6



    Mg_fit=line_fitter(sp, "MgII", i, guesses_MgII, limits_MgII, limited_MgII, tied_MgII, xmin_MgII, xmax_MgII,offset=-0.05, magorder=mag_order)
    print "FWHM_Mg=", sp.specfit.measure_approximate_fwhm()
    FWHM_Mg=sp.specfit.measure_approximate_fwhm()
    FWHMV_Mg=wl_to_kms(FWHM_Mg,2800)
    smooth_f=FWHMV_Mg/600.0 
    sigmav_Mg=FWHMV_Mg/(2.355)
    
    sigma=kms_to_wl(np.sqrt(sigmav_Mg*sigmav_Mg-600.0*600.0),3800) 
    balmer_tot=gauss_conv(balmer_tot, sigma, order=0, output=None, mode='reflect', cval=0.0)
    print "smooth", smooth_f,sigmav_Mg

    balmer_lines.data=smooth(balmer_lines.data,smooth_f , downsample=False)
    #balmer_lines.smooth(smooth_f)
    #balmer_tot=smooth(balmer_tot,smooth_f , downsample=False)
    
    #index=np.argmin( np.abs( lambda_balmer - sp.xarr) )
    #index_hi=np.argmin( np.abs( 4000.0 - sp.xarr) )
    balmer_data=balmer_tot
    #np.zeros(len(balmer_template.data))
    #balmer_data[:index]=balmer_template.data[:index]
    #balmer_data[index:index_hi]=balmer_lines.data[index:index_hi]
    xmin_bb=3640
    xmax_bb=3650
    

    pylab.figure()
    pylab.plot(backup.xarr,backup.data)
    pylab.plot(balmer_template.xarr[:index],balmer_template.data[:index],balmer_lines.xarr[index:],balmer_lines.data[index:] )
    pylab.plot(sp.xarr,sp.data)
    pylab.savefig( "./plots/balmer_fit_"+str(i)+".pdf" )
    
    



    
    #fe_fitter("C",i,xmin_CIII,xmax_CIII)
    Si_fit=line_fitter(sp,"SiIV", i,guesses_SiIV, limits_SiIV, limited_SiIV, tied_SiIV, xmin_SiIV, xmax_SiIV, magorder=mag_order)
    print "FWHM_Si=", sp.specfit.measure_approximate_fwhm()
    C_fit=Si_fit*0.0
    
    C_fit=line_fitter(sp, "C", i,guesses_C, limits_C, limited_C, tied_C, xmin_C, xmax_C,magorder=mag_order)
    
    print "FWHM_C=", sp.specfit.measure_approximate_fwhm()


    
    
    xmin_fe=4600
    xmax_fe=4800
    template_file="./templates/fe2_Op2.dat"
    fe_template = pyspeckit.Spectrum(template_file)
    fe_template.xarr.units='angstroms'
    fe_template.xarr.xtype='wavelength'
    lambda_fe=4500
    index_sp=np.argmin( np.abs( lambda_fe - sp.xarr) )
    index_fe=np.argmin( np.abs( lambda_fe - fe_template.xarr) )
    flux=np.mean(sp.data[ index_sp - 15 : index_sp + 15 ])
    fe_template.data*=np.abs(flux/fe_template.data[index_fe])
    pylab.rcParams["figure.figsize"]=16,6
    #fe_Hb=blue_bump_fitter(sp,"optical",i,xmin_fe,xmax_fe,fe_template,mag_order)
    pylab.rcParams["figure.figsize"]=8,6
    #print limits_Hbeta
    Hbeta_fit=line_fitter(sp, "Hbeta", i, guesses_Hbeta, limits_Hbeta, limited_Hbeta, tied_Hbeta, xmin_Hbeta, xmax_Hbeta,offset=-0.05, magorder=mag_order)
    print "FWHM_Hbeta=", sp.specfit.measure_approximate_fwhm()
    Halpha_fit=line_fitter(sp, "Halpha", i, guesses_Halpha, limits_Halpha, limited_Halpha, tied_Halpha, xmin_Halpha, xmax_Halpha,offset=-0.05, magorder=mag_order)
    print "FWHM_Halpha=", sp.specfit.measure_approximate_fwhm()
    


    
    print len(Halpha_fit)
    print len(Hbeta_fit)
    print len(C_fit)
    print len(Mg_fit)
    print len(fe_fit)
    print len(balmer_template.data)
    total=continuous +  Halpha_fit + Hbeta_fit + C_fit + fe_fit  + balmer_data + Si_fit + Mg_fit

    pylab.rcParams["figure.figsize"]=16,6


    pylab.figure()
    pylab.ylim(ymin=-1.3*np.abs(copy1.data.min()),ymax=1.1*copy1.data.max())
    pylab.xlim(xmin=wlmin,xmax=wlmax)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.yscale('linear')
    pylab.plot(sp.xarr,copy3.data,'k')
    pylab.plot(sp.xarr, total,'r')
    total=continuous + fe_fit  + balmer_data 
    pylab.plot(sp.xarr, total,'b')
    total=continuous + balmer_data 
    pylab.plot(sp.xarr, total,'g')
    total=continuous
    pylab.plot(sp.xarr, total,'gray')
    
    plot_file="./plots/superposition_" + str(i) + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)
    
    
    pylab.figure()
    xmin=2000
    xmax=4000
    index_min=np.argmin(np.abs(xmin-sp.xarr))
    index_max=np.argmin(np.abs(xmax-sp.xarr))
    pylab.ylim(ymax=1.1*sp.data[index_max])

    pylab.xlim(xmin=xmin,xmax=xmax)
    pylab.ylabel(r'$10^{'+str(mag_order)+'}$ erg s$^{-1}$  $cm^{-1}$')
    pylab.xlabel(r'$\AA$')
    pylab.yscale('linear')
    pylab.plot(sp.xarr,copy3.data,'k')
    pylab.plot(sp.xarr, total,'r')
    total=continuous + fe_fit  + balmer_tot
    pylab.plot(sp.xarr, total,'b')
    #    total=continuous + fe_fit + balmer_tot1
    #   pylab.plot(sp.xarr, total, "pink")
    
    total=continuous + balmer_data 
    pylab.plot(sp.xarr, total,'g')

    total=continuous
    pylab.plot(sp.xarr, total,'gray')
    
    plot_file="./plots/mg+fe_" + str(i) + ".png"
    #sp.plotter.figure.savefig(plot_file)
    pylab.savefig(plot_file)

    pylab.rcParams["figure.figsize"]=8,6


    print L_cont







