#!/usr/bin/python



# -*- coding: utf-8 -*-

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa
#Standalone_Taufit_simu.py

import argparse
import re
import os, sys
import pypsr_standalone as psr
import matplotlib.pyplot as plt
import lmfit
from lmfit import Model, conf_interval, printfuncs
from lmfit import minimize, Parameter, Parameters, fit_report
from lmfit.models import LinearModel, PowerLawModel, ExponentialModel, QuadraticModel
import numpy as np
from scipy import special
from scipy import stats
import DataReadIn as dri

from matplotlib import rc
rc('text', usetex=True)

#"""Read and print header information"""
"""Define options to the script"""
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename',
                    help="Provide the pathway to the data files")
parser.add_argument('-m','--method',
                    help="Choosing method to fit data or simulation. Choose between 'onedim', 'iso', 'aniso','postfold', 'isoplusonedim'")                   
parser.add_argument('-dc','--datacycle',
                    help="The type of data. Choose from comm., census, cycle5. Only used in creating filenames.")


args = parser.parse_args()

"""Allocate variable names to the parsed options"""
filepath = args.filename
meth = args.method
datac = args.datacycle

"""Create folder to save to"""
newpath = r'./SummaryPlots'
if not os.path.exists(newpath):
    os.makedirs(newpath)


print "\n Reading in data \n"
pulsars = []
pulsar, nch, nbins,nsub, lm_rms, tsub = dri.read_headerfull(filepath)
pulsars = [pulsar for i in range(nsub)]
print3 = "RMS: %f" %lm_rms
 
print0 = "Pulsar name: %s" %pulsar
print1 = "Number of channels: %d" %nch
print2 = "Number of bins: %d" %nbins

for k in range(4):
    print eval('print{0}'.format(k))
    

profilexaxis = np.linspace(0,tsub,nbins)

#for i in range(15):
#    plt.close(i)


## Create subplot structure
sp = 12 #number of subplots per figure
#plt.figure(figsize=(16,10))  

#Fit profile using tau_fitter

obtainedtaus = []
lmfittausstds = []

bestparamsall = []
bestparams_stdall = []
redchis = []

freqmsMHz =[]
noiselessmodels =[]
results = []
comp_rmss = []
comp_fluxes= []
comp_SNRs =[]
datas = []



#trainlength = 4

halfway = nbins/2.


for i in range(nsub):
    data, freqm = dri.read_data(filepath,i,nbins)
    freqmsMHz.append(freqm) 
    peakbin = np.argmax(data)
    shift = int(halfway)-int(peakbin)
    data = np.roll(data,shift)
    print "I'm rolling the data to ensure peak is in the middle"
    comp_rms = psr.find_rms(data,nbins)
   
    if meth is None:
		print "No fitting method was chosen. Will default to an isotropic fitting model. \n Use option -m with 'onedim' or 'aniso' to change."
		result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_fitter(data,nbins)
    elif meth in ('iso','Iso','ISO'):
		result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_fitter(data,nbins)
    elif meth in ('onedim','1D','Onedim'):
		result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_1D_fitter(data,nbins)
    elif meth in ('aniso','Aniso','ANISO'):
		noiselessmodel, besttau, taustd, besttau2, taustd2, bestparams, bestparams_std, redchi = psr.tau_ani_fitter(data,nbins)
    elif meth in ('multi'):
              if i == 0:
#        		    bw1 = raw_input("Provide bin-window for 1st the peak (e.g: [20,30]): ")
#        		    bw1 = eval(bw1)
#        		    bw2 = raw_input("Provide bin-window start and end 2nd peak (e.g: [35,45]): ")
#        		    bw2 = eval(bw2)
        		    bw1, bw2 = np.linspace(660,720,2), np.linspace(725,750,2)        
        		    peak1 = np.max(data[bw1[0]: bw1[-1]])
        		    peak2 = np.max(data[bw2[0]: bw2[-1]])
        		    maxmean = 2.0*bw2[-1]
              else:
                     bw1,bw2 = bw1,bw2
                     peak1, peak2, maxmean = peak1, peak2, maxmean
              nms, bt, tstd, bp, bpstd, chis = [], [], [], [], [], []
              for k in range(0,len(bw1)):
                  for j in range(0,len(bw2)):
        			noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_fittermulti(data,nbins,bw1[k],bw2[j], peak1, peak2, maxmean)
        			nms.append(noiselessmodel)
        			bt.append(besttau)
        			tstd.append(taustd)
        			bp.append(bestparams)
        			bpstd.append(bestparams_std)
        			chis.append(redchi)
        			print j*k
              tstd = np.array(tstd)
              nonz = np.nonzero(tstd)
              nonz_min = np.argmin(tstd[np.nonzero(tstd)])
              tstd_ind = np.array(nonz).flatten()[nonz_min]     
              noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = nms[tstd_ind], bt[tstd_ind], tstd[tstd_ind], bp[tstd_ind], bpstd[tstd_ind], chis[tstd_ind]               
    else:
             print "Incorrect fitting method. Choose from iso, onedim, aniso, postfold"    


    """ INPUT SECTION ENDS """

            
    comp_flux = psr.find_modelflux(noiselessmodel,nbins)
    comp_SNR_model = psr.find_peaksnr(noiselessmodel,comp_rms)

    print 'SNR (from model): %.2f' % comp_SNR_model
    comp_SNR =  psr.find_peaksnr_smooth(data,comp_rms)
    print 'SNR (from data): %.2f' % comp_SNR

    obtainedtaus.append(besttau)    
    lmfittausstds.append(taustd)
    bestparamsall.append(bestparams)
    bestparams_stdall.append(bestparams_std)
    redchis.append(redchi)
    
        
    noiselessmodels.append(noiselessmodel)
    results.append(result)
    comp_rmss.append(comp_rms)
    comp_fluxes.append(comp_flux)
    comp_SNRs.append(comp_SNR)
    datas.append(data)


for k in range(0,4):
    print eval('print{0}'.format(k))
            

data_highsnr = datas[0]
freqms_highsnr = np.array(freqmsMHz[0])/1000.
freqMHz_highsnr = freqmsMHz[0]
taus_highsnr = obtainedtaus[0]
lmfitstds_highsnr = np.array(lmfittausstds)[0]
model_highsnr = noiselessmodels[0]
fluxes_highsnr = comp_fluxes[0]
rms_highsnr = comp_rmss[0]
redchis_highsnr = redchis[0]


"""Array with all the other fitting parameters: sigma, A, etc."""
bestpT = np.transpose(bestparamsall)
bestpT_std = np.transpose(bestparams_stdall)

bestpT_highSNR = bestpT
bestpT_std_highSNR = bestpT_std

print6 = "Number of plots %d" %(nsub)
npch = nsub

print eval('print{0}'.format(6))

tbs = tsub/nbins
    
"""Plotting starts"""

plt.close('all')

alfval = 1.0
alfval2 = 0.2
markr = 'k*'
prof = 'r-'
textpos = 0.8
textpos2 = 5
    
##PLOT PROFILES##

numFig = 1

profiles = []
subplotcount = 1
figg = plt.figure(numFig,figsize=(16,10))
figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05)    
plt.subplot(2,2,subplotcount)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(profilexaxis,data_highsnr,'y',alpha = 0.25)
plt.plot(profilexaxis,model_highsnr,prof, alpha = 0.7)
plt.title('%s at %.1f MHz, Tsub: %.2f sec' %(pulsar, freqMHz_highsnr, tsub))
plt.text(0.6*tsub,0.8*np.max(data_highsnr), r'$\tau: %.2f \pm %.2f$ ms' %(taus_highsnr*tbs*1000, lmfitstds_highsnr*tbs*1000),fontsize=14)
plt.ylim(ymax=1.3*np.max(data_highsnr))
plt.xlim(xmax=tsub)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('time (s)',fontsize=14)
plt.ylabel('normalized intensity',fontsize=14)
        

"""Compute the KS Test D value and probability of the residuals following a normal distribution"""

resdata = data_highsnr - model_highsnr
resnormed = (resdata-resdata.mean())/resdata.std()
KSd, KSp = stats.kstest(resnormed, 'norm')

#print "Probability of residuals being Gaussian: %.4f" % KSp
#print "Reduced Chi square: %.2f" %  np.mean(redchis_highsnr/np.power(rms_highsnr,2))

##PLOT RESIDUALS

subplotcount += 1

plt.subplot(2,2,subplotcount)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(profilexaxis,resdata,'b',alpha = 0.25)
plt.title('%s at %.1f MHz' %(pulsar, freqMHz_highsnr))
plt.xlim(xmax=tsub)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('time (s)',fontsize=14)
plt.ylabel('residuals',fontsize=14)

##PLOT RESIDUALS HISTOGRAM

subplotcount += 1

plt.subplot(2,2,subplotcount)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.hist(resdata,facecolor='b',bins=10)
plt.title('%s at %.1f MHz' %(pulsar, freqMHz_highsnr))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('time (s)',fontsize=14)
plt.ylabel('counts',fontsize=14)


###PLOT FLUX##
#
xaxis = np.linspace(1,npch,npch)
#
#plt.figure(numfig, figsize=(16,10))
#plt.subplot(2,2,2)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.plot(xaxis,fluxes_highsnr,markr, alpha = alfval)
#plt.title('%s' %(pulsar))
#plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)
#plt.xlabel('Tsubs (count)',fontsize=16)
#plt.ylabel(r'Calibrated flux (mJy)',fontsize=16)


lmfittausstds = np.array(lmfittausstds)
obtainedtaus = np.array(obtainedtaus)


obtainedtausec = obtainedtaus*tbs
lmfitstdssec_highsnr = lmfitstds_highsnr*tbs

freqms = np.array(freqmsMHz)/1000. #(Use this one to compare the alpha outcome of all the freq channels with only the high SNR/low tau error ones.)

"""The associated high SNR arrays are"""
taus_highsnr = np.array(taus_highsnr)
taussec_highsnr = taus_highsnr*tbs


for i in range(nsub):
    print'Tau (ms): %.2f' %(1000*taussec_highsnr)

print9 = 'Tsub = %.2f sec' %tsub
print eval('print{0}'.format(9))
    
##PLOT TAU##  

#plt.figure(numFig, figsize=(16,10))
#plt.subplot(2,2,3)
##plt.errorbar(xaxis,taussec_highsnr,yerr=lmfitstds_highsnr*tbs,fmt=markr,markersize=9.0,capthick=2,linewidth=1.5, alpha = alfval)
#plt.plot(xaxis,taussec_highsnr,markr,markersize=9.0,linewidth=1.5, alpha = alfval)
##plt.plot(1000.*freqms_highsnr, powouttau_highsnr.best_fit, 'k-', alpha=alfval)
#plt.plot(xaxis,fluxes_highsnr,markr, alpha = alfval)
#plt.title('%s' %(pulsar))
#plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)
#plt.xlabel('Tsubs (count)',fontsize=16)
#plt.ylabel(r'$\tau$ (sec)',fontsize=16)



### PLOT CHI ##  

#plt.figure(numFig, figsize=(16,10))
#plt.subplot(2,2,4)
#plt.plot(xaxis, redchis_highsnr/np.power(rms_highsnr,2), markr,alpha=alfval,markersize = 12)
#plt.title(r'Reduced $\chi^2$ values')
#plt.yticks(fontsize=10)
#plt.xticks(fontsize=12)
#plt.xlabel('Tsubs (count)',fontsize=16)
#plt.ylabel(r'$\chi^2$',fontsize=16)



##PLOT SIGMA##  #

#figg = plt.figure(numFig, figsize=(16,10))
#figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05)    
#plt.subplot(2,2,4)
#plt.plot(xaxis,bestpT_highSNR[0]*tbs,'m*',markersize=9.0,linewidth=1.5,alpha=alfval)
##plt.errorbar(xaxis,bestpT_highSNR[0]*tbs,yerr = bestpT_std_highSNR[0]*tbs, fmt = 'm*',markersize=9.0,capthick=2,linewidth=1.5,alpha=alfval)
#plt.ylabel(r'$\sigma$ (sec)')
#plt.yticks(fontsize=11)
#plt.xticks(fontsize=10)
#plt.xlabel('Tsubs (count)',fontsize=16)
#plt.legend(fontsize = 9, loc='best')
#
#
### PLOT A ##
#
plt.figure(numFig,figsize=(16,10))
plt.subplot(2,2,4)
#plt.errorbar(xaxis, bestpT_highSNR[2], yerr = bestpT_std_highSNR[2], fmt = 'g*',markersize=9.0,capthick=2,linewidth=1.5, alpha=alfval)
plt.plot(xaxis, bestpT_highSNR[2],'g*',markersize=9.0,linewidth=1.5, alpha=alfval)
plt.title('Amplitude')
plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xlabel('Tsubs (count)',fontsize=16)
plt.ylabel(r'$A$',fontsize=16)


print 'Sigma: %.2f ms' % (bestpT_highSNR[0]*tbs*1000)
print 'Amp: %.2f ' % (bestpT_highSNR[2])

savename = re.split("Beam2_",filepath)[-1]
savename = savename[0:-9]

Summaryplot = '%s_%s.png'  % (savename,meth)
picpathtau = newpath
fileoutputtau = os.path.join(picpathtau,Summaryplot)
plt.savefig(fileoutputtau, dpi=150)
print 'Saved %s in %s' %(Summaryplot,newpath)
#plt.close()
