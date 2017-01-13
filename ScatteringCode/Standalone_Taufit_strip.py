#!/usr/bin/python



# -*- coding: utf-8 -*-

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa
#Standalone_Taufit_simu.py

import argparse
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
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#"""Read and print header information"""
"""Define options to the script"""
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename',
                    help="Provide the pathway to the data files")
parser.add_argument('-p','--period',type=float,
                    help="Provide the period of the pulsar in seconds")
parser.add_argument('-c','--counts',type=int,
                    help="Provide the number of noise realisations (interger)")
parser.add_argument('-m','--method',
                    help="Choosing method to fit data or simulation. Choose between 'onedim', 'iso', 'aniso','postfold', 'isoplusonedim'")                   
parser.add_argument('-dc','--datacycle',
                    help="The type of data. Choose from comm., census, cycle5. Only used in creating filenames.")


args = parser.parse_args()

"""Allocate variable names to the parsed options"""
filepath = args.filename
pulseperiod = args.period
count = args.counts
meth = args.method
datac = args.datacycle

"""Create folder to save to"""
newpath = r'./SummaryPlots'
if not os.path.exists(newpath):
    os.makedirs(newpath)


print "\n Reading in data \n"
pulsars = []
pulsar, nch, nbins,nsub, lm_rms = dri.read_headerfull(filepath)
pulsars = [pulsar for i in range(nsub)]
print3 = "RMS: %f" %lm_rms
 
print0 = "Pulsar name: %s" %pulsar
print1 = "Number of channels: %d" %nch
print2 = "Number of bins: %d" %nbins

for k in range(4):
    print eval('print{0}'.format(k))
    

#Find pulseperiod from list if it wasn't parsed

pulsarBnamelist = ['B0037+56','B0114+58','B0540+23','B0611+22','B0740-28','B1848+12','B1907+10','B1911-04','B1915+13','B1920+21','B1933+16','B2255+58','B2303+30']
pulsarnamelist = ['J0040+5716','J0117+5914','J0543+2329','J0614+2229', 'J0742-2822','J1851+1259','J1909+1102','J1913-0440','J1917+1353','J1922+2110','J1935+1616','J2257+5909','J2305+3100']
pulsarperiodlist = [1.118225, 0.101439, 0.245975, 0.33496, 0.166762, 1.205303, 0.283641, 0.825936, 0.194631, 1.077924, 0.358738, 0.368246, 1.575886]


if pulseperiod == None: 
   	pulsarN = pulsar[0:8]
   	print pulsarN
   	pulseind = pulsarBnamelist.index(pulsarN)
   	pulseperiod = pulsarperiodlist[pulseind]
else:
   pulseperiod = pulseperiod


profilexaxis = np.linspace(0,pulseperiod,nbins)

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
            

data_highsnr = datas
freqms_highsnr = np.array(freqmsMHz)/1000.
freqMHz_highsnr = freqmsMHz
taus_highsnr = obtainedtaus
lmfitstds_highsnr = np.array(lmfittausstds)
model_highsnr = noiselessmodels
fluxes_highsnr = comp_fluxes
rms_highsnr = comp_rmss
redchis_highsnr = redchis


"""Array with all the other fitting parameters: sigma, A, etc."""
bestpT = np.transpose(bestparamsall)
bestpT_std = np.transpose(bestparams_stdall)

bestpT_highSNR = bestpT
bestpT_std_highSNR = bestpT_std

print6 = "Number of plots %d" %(nsub)
npch = nsub

print eval('print{0}'.format(6))

"""Calculate fits for parameters sigma and mu"""
"""Shift data based on mu-trend (DM trend?)"""

pbs = pulseperiod/nbins
    
"""Plotting starts"""

plt.close('all')
alfval = 1.0
alfval2 = 0.2
markr = 'k*'
prof = 'r-'
textpos = 0.8
textpos2 = 5
    
##PLOT PROFILES##

totFig = (nsub+5)/sp + 1


profiles = []
for j in range(nsub):
    if j+1 == sp:
        numFig = (j+1)/sp
        subplotcount = sp
    else: 
        numFig = (j+1)/sp + 1
        if j+1 < sp:
            subplotcount = j+1
        else: 
            subplotcount = j+1 - sp
    figg = plt.figure(numFig,figsize=(16,10))
    figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05)    
    plt.subplot(3,4,subplotcount)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(profilexaxis,data_highsnr[j],'y',alpha = 0.25)
    plt.plot(profilexaxis,model_highsnr[j],prof, alpha = 0.7)
    plt.title('%s at %.1f MHz' %(pulsars[j], freqMHz_highsnr[j]))
    plt.annotate(r'$\tau: %.4f \pm %.2e$ sec' %(taus_highsnr[j]*pulseperiod/nbins, lmfitstds_highsnr[j]*pulseperiod/nbins),xy=(np.max(profilexaxis),np.max(data_highsnr[j])),xycoords='data',xytext=(0.4,textpos),textcoords='axes fraction',fontsize=12)
    plt.ylim(ymax=1.1*np.max(data_highsnr[j]))
    plt.xlim(xmax=pulseperiod)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('time (s)')
    plt.ylabel('normalized intensity')
        

if npch >= sp:
    subpl_ind = int(npch-sp) +1
else:
    subpl_ind = int(npch) + 1

if npch + 1 == sp:
    numfig = (npch+1)/sp
else:
    numfig = (npch+1)/sp + 1

##PLOT FLUX##

xaxis = np.linspace(1,npch,npch)

plt.figure(numfig, figsize=(16,10))
plt.subplot(3,4,subpl_ind)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(xaxis,fluxes_highsnr,markr, alpha = alfval)
plt.title('%s' %(pulsar))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Tsubs (count)',fontsize=16)
plt.ylabel(r'Calibrated flux (mJy)',fontsize=16)

lmfittausstds = np.array(lmfittausstds)
obtainedtaus = np.array(obtainedtaus)


obtainedtausec = obtainedtaus*pulseperiod/nbins
lmfitstdssec_highsnr = lmfitstds_highsnr*pulseperiod/nbins

freqms = np.array(freqmsMHz)/1000. #(Use this one to compare the alpha outcome of all the freq channels with only the high SNR/low tau error ones.)

"""The associated high SNR arrays are"""
taus_highsnr = np.array(taus_highsnr)
taussec_highsnr = taus_highsnr*pulseperiod/nbins


print9 = 'pulseperiod = %.6f' %pulseperiod
for i in range(nsub):
    print'%d Tau (ms): %.2f' %(i, 1000*taussec_highsnr[i])

print eval('print{0}'.format(9))
    
##PLOT TAU##  
    
if subpl_ind >= sp:
    subpl_ind2 = int(subpl_ind-sp) +1
else:
    subpl_ind2 = int(subpl_ind) + 1

if npch + 2 == sp:    
    numfig2 = (npch+2)/sp
else:
    numfig2 = (npch+2)/sp + 1

plt.figure(numfig2, figsize=(16,10))
plt.subplot(3,4,subpl_ind2)
#plt.errorbar(xaxis,taussec_highsnr,yerr=lmfitstds_highsnr*pulseperiod/nbins,fmt=markr,markersize=9.0,capthick=2,linewidth=1.5, alpha = alfval)
plt.plot(xaxis,taussec_highsnr,markr,markersize=9.0,linewidth=1.5, alpha = alfval)
#plt.plot(1000.*freqms_highsnr, powouttau_highsnr.best_fit, 'k-', alpha=alfval)
plt.plot(xaxis,fluxes_highsnr,markr, alpha = alfval)
plt.title('%s' %(pulsar))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Tsubs (count)',fontsize=16)
plt.ylabel(r'$\tau$ (sec)',fontsize=16)



### PLOT CHI ##  

if subpl_ind2 >= sp:
    subpl_ind3 = int(subpl_ind2-sp) +1
else:
    subpl_ind3 = int(subpl_ind2) + 1
    
if npch + 3 == sp:    
    numfig3 = (npch+3)/sp
else:
    numfig3 = (npch+3)/sp + 1

plt.figure(numfig3, figsize=(16,10))
plt.subplot(3,4,subpl_ind3)
plt.plot(xaxis, redchis_highsnr/np.power(rms_highsnr,2), markr,alpha=alfval,markersize = 12)
plt.title(r'Reduced $\chi^2$ values')
plt.yticks(fontsize=10)
plt.xticks(fontsize=12)
plt.xlabel('Tsubs (count)',fontsize=16)
plt.ylabel(r'$\chi^2$',fontsize=16)

"""Compute the KS Test D value and probability of the residuals following a normal distribution"""
KSs = np.zeros((npch,2))
ADs = np.zeros((npch,2))
for i in range(npch):
    resdata = data_highsnr[i] - model_highsnr[i]
    resnormed = (resdata-resdata.mean())/resdata.std()
    KSd, KSp = stats.kstest(resnormed, 'norm')
    KSs[i,0] = KSd
    KSs[i,1]= KSp
    
#    aa,bb,cc = stats.anderson(resnormed, 'norm')
#    print aa,bb,cc

resmeanP = np.mean(KSs[:,1])
print "Mean probability of residuals being Gaussian: %.4f" % resmeanP

print "Mean reduced Chi square: %.2f" %  np.mean(redchis_highsnr/np.power(rms_highsnr,2))




##PLOT SIGMA##  #

if subpl_ind3 >= sp:
    subpl_ind4 = int(subpl_ind3-sp)+1
else:
    subpl_ind4 = int(subpl_ind3)+1

if npch + 4 == sp:    
    numfig4 = (npch+4)/sp
else:
    numfig4 = (npch+4)/sp + 1


figg = plt.figure(numfig4, figsize=(16,10))
figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05)    
plt.subplot(3,4,subpl_ind4)
plt.plot(xaxis,bestpT_highSNR[0]*pbs,'m*',markersize=9.0,linewidth=1.5,alpha=alfval)
#plt.errorbar(xaxis,bestpT_highSNR[0]*pbs,yerr = bestpT_std_highSNR[0]*pbs, fmt = 'm*',markersize=9.0,capthick=2,linewidth=1.5,alpha=alfval)
plt.ylabel(r'$\sigma$ (sec)')
plt.yticks(fontsize=11)
plt.xticks(fontsize=10)
plt.xlabel('Tsubs (count)',fontsize=16)
plt.legend(fontsize = 9, loc='best')


## PLOT A ##

if subpl_ind4 >= sp:
    subpl_ind5 = int(subpl_ind4-sp)+1
else:
    subpl_ind5 = int(subpl_ind4)+1

if npch + 5 == sp:    
    numfig5 = (npch+5)/sp
else:
    numfig5 = (npch+5)/sp + 1


plt.figure(numfig5,figsize=(16,10))
plt.subplot(3,4,subpl_ind5)
#plt.errorbar(xaxis, bestpT_highSNR[2], yerr = bestpT_std_highSNR[2], fmt = 'g*',markersize=9.0,capthick=2,linewidth=1.5, alpha=alfval)
plt.plot(xaxis, bestpT_highSNR[2],'g*',markersize=9.0,linewidth=1.5, alpha=alfval)
plt.title('Amplitude')
plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xlabel('Tsubs (count)',fontsize=16)
plt.ylabel(r'$A$',fontsize=16)


for i in range(numfig5):
    k = numfig5 - i ##reverse the order
    Summaryplot = '%s_%s_%s_%d.png'  % (pulsar,datac,meth,k)
    picpathtau = newpath
    fileoutputtau = os.path.join(picpathtau,Summaryplot)
    plt.savefig(fileoutputtau, dpi=150)
    print 'Saved %s in %s' %(Summaryplot,newpath)
    plt.close()
