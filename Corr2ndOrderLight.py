import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import glob
import pickle

ispec_dir = '/Users/chimamcgruder/Research/'
sys.path.append(ispec_dir+'TwinPlanets/')
sys.path.append(ispec_dir+'iSpec_v20201001')
import ispec 
import AnalyzeDat as AD  

#### PATH info #### 
BasePath = '/Users/chimamcgruder/Research/HydraDump/WASP96/'
BasePath1, BasePath2 = BasePath+'ut170804_05/',BasePath+'ut171108_09/'
tickFont, Font = 12, 18
def DegradeRes(wave, flux, initiRes, WaveRange=None, finalRes=10, filter_regions=None): #To degrade the resolution of the spectrum, to determine the continuum of the spectrum
    # filter_regions = regions that should be ignorned when making the smoothed spectra. Could be a list of list or a .txt file
    if WaveRange:
        WavRng = np.where((WaveRange[0]<wave) & (wave<WaveRange[1]))
        wave, flux = wave[WavRng], flux[WavRng]
    if filter_regions:
        if type(filter_regions) == str: # then it's a file
            DnWav, UpWav = np.loadtxt(filter_regions, unpack=True)
        if type(filter_regions) == list: #then N x 2 list
            y = np.array(filter_regions) #convert to 2 x N array so in the same structure as np.loadtxt
            z = np.transpose(y)#convert to 2 x N array so in the same structure as np.loadtxt
            DnWav, UpWav = z[0], z[1]
        for fr in range(len(UpWav)):
            BadFrames = np.where((DnWav[fr]<wave) & (wave<UpWav[fr]))
            wave, flux = np.delete(wave,BadFrames), np.delete(flux,BadFrames)
    col1 = fits.Column(name='WAVE', format='D', array=wave) #original wavelength grid
    col2 = fits.Column(name='Flux', format='D', array=flux)
    cols = fits.ColDefs([col1, col2])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr['TTYPE1'], prihdr['TTYPE2'] = 'WAVE',"FLUX"
    prihdr['TFORM1'], prihdr['TFORM2'] = '313111D',"313111D"
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    TempFilName = '/TempSpec_Test.fits'
    New_res_path='/Users/chimamcgruder/Research/TwinPlanets/WASP96/IMACS'
    thdulist.writeto(New_res_path+TempFilName, overwrite=True) #to tempirarily save the spectra in a fits file, in the right formate for iSpec to read
    thdulist.close()
    ResDeSpec, OGspec = AD.degrade_resolution(New_res_path+TempFilName, from_res=initiRes, to_res=finalRes)
    os.remove(New_res_path+TempFilName)
    return ResDeSpec, OGspec

#load GG495 throughput
Wave, through = np.loadtxt(BasePath+'GG495filter.txt', unpack=True)
Wave, through = Wave*10, through/100.0 #to convert from nm to angstroms #and % transmission to fractional transmission

plt.figure(figsize=(18, 4))
########### load spectra and okit spectra 
y1 = fits.getdata(BasePath1+'data_reductions17_120_noOpt/final_spec_WASP-96.fits')
wave1, flux1 = y1[0,:], np.mean(y1[1:,:], axis=0)
Chip2idx = np.where(wave1>8850.0)[0]
flux1[Chip2idx] = flux1[Chip2idx]*.85 #difference in chip gaps was emperically found to be about .85
ThroInterp1 = np.interp(wave1, Wave, through)
FilteredFlux1 = flux1*ThroInterp1
Mx1 = np.argmax(ThroInterp1) #index location of max throughput of filter. Use this location to maximize normalize function

y2 = fits.getdata(BasePath2+'data_reductions8_95_noOpt/final_spec_WASP-96.fits')
wave2, flux2 = y2[0,:], np.mean(y2[1:,:], axis=0)
Chip2idx = np.where(wave2>8850.0)[0]
flux2[Chip2idx] = flux2[Chip2idx]*.85 #difference in chip gaps was emperically found to be about .85
ThroInterp2 = np.interp(wave2, Wave, through)
Mx2 = np.argmax(ThroInterp2) #index location of max throughput of filter. Use this location to maximize normalize function

plt.subplot2grid((1,8), (0,0), colspan=3)
normFlux1, normFiltFlux1, normFlux2 = flux1/flux1[Mx1], FilteredFlux1/FilteredFlux1[Mx1], flux2/flux2[Mx2]
Gap1, Gap2 = np.where((8788.8 < wave1) & (wave1 < 8892.3))[0], np.where((8788.8 < wave2) & (wave2 < 8892.3))[0] 
normFlux1[Gap1], normFiltFlux1[Gap1], normFlux2[Gap2] = np.nan, np.nan, np.nan
plt.plot(wave1/10, normFlux1, 'g', label = 'ob1 (UT170804)')
plt.plot(wave2/10, normFlux2, 'b', label = 'ob2 (UT171108)')
plt.plot(wave1/10, normFiltFlux1, 'r', label = 'ob1*GG495')
plt.ylabel("Normalized flux", fontweight='bold', fontsize=Font)
plt.xlabel("Wavelength [nm]", fontweight='bold', fontsize=Font)
plt.xticks(fontweight='bold', fontsize=tickFont)
plt.yticks(fontweight='bold', fontsize=tickFont)
plt.legend()
plt.xlim((490,940))


########### to fit the light curve's contiumm and plot it
wav_rang = [5000.0, 9400.0] #wavelegth range where the continuum fit is valid

#For some reason in the .py file the continua a plotted weird. Don't deal with it and just load what the notebook did

plt.subplot2grid((1,8), (0,5), colspan=3)
# Observation 1
ResDeSpec, OGspec = DegradeRes(wave1, flux1, 700, finalRes=10, WaveRange=[4900.0, 9900.0], filter_regions = 'W96_ut170804_screensW96.txt')
fluxNormilization = flux1[np.where(wave1==5683.75)[0]] #normalization term to set both observations to have the same blue flux
wave_used, flux_used = OGspec['waveobs']*10, OGspec['flux']
NormFLUX1 = flux1/fluxNormilization
Ratio = np.max(ResDeSpec['flux'])/np.max(flux1) #to get ratio of maxes of the continuum and spectra so the normalized continuum is mapped on top of the normalized spectra 
ResWav1, ResNormFlux1 = ResDeSpec['waveobs']*10, ResDeSpec['flux']/(fluxNormilization*Ratio)#[np.where(wave== 5683.75)[0]]
WavRng = np.where((wav_rang[0]<ResWav1) & (ResWav1<wav_rang[1]))
cResWav, cResNormFlux = np.copy(ResWav1), np.copy(ResNormFlux1)
cutResWav1, cutResNormFlux1 = cResWav[WavRng], cResNormFlux[WavRng]
# plt.plot(cutResWav1, cutResNormFlux1, 'r--')
X = np.load('continum_UT170804.npy')
Xx = np.load('UT170804_spectrum.npy')
plt.plot(X[0,:]/10, X[1,:], 'g--', label='ob1 (UT170804)')
Gap = np.where((8788.8 < Xx[0,:]) & (Xx[0,:] < 8892.3))[0]
Xx[1,:][Gap] = np.nan
plt.plot(Xx[0,:]/10, Xx[1,:], 'g.', markersize=1)
# plt.plot(wave1, NormFLUX1, 'g.', markersize=1, label='observation UT170804')

# Observation 2
ResDeSpec, OGspec = DegradeRes(wave2, flux2, 700, finalRes=10, WaveRange=[4900.0, 9900.0], filter_regions = 'W96_ut171108_screensW96.txt')
fluxNormilization = flux2[np.where(wave2== 5683.75)[0]] #normalization term to set both observations to have the same blue flux
wave_used, flux_used = OGspec['waveobs']*10, OGspec['flux']
NormFLUX2 = flux2/fluxNormilization
Ratio = np.max(ResDeSpec['flux'])/np.max(flux2) #to get ratio of maxes of the continuum and spectra so the normalized continuum is mapped on top of the normalized spectra 
ResWav, ResNormFlux = ResDeSpec['waveobs']*10, ResDeSpec['flux']/(fluxNormilization*Ratio)
WavRng = np.where((wav_rang[0]<ResWav) & (ResWav<wav_rang[1]))
cResWav, cResNormFlux1= np.copy(ResWav), np.copy(ResNormFlux)
cutResWav2, cutResNormFlux2 = cResWav[WavRng], cResNormFlux[WavRng]
ReNormTerm = (cutResNormFlux1[0]/cutResNormFlux2[0])
# plt.plot(cutResWav2, cutResNormFlux2*ReNormTerm, 'b')
Y = np.load('UT171108_spectrum.npy')
Yy = np.load('continum_UT171108.npy')
plt.plot(Yy[0,:]/10, Yy[1,:], 'b', label='ob2 (UT171108)')
Gap = np.where((8788.8 < Y[0,:]) & (Y[0,:] < 8892.3))[0]
Y[1,:][Gap] = np.nan
plt.plot(Y[0,:]/10, Y[1,:], 'b.', zorder=0, markersize=1)
plt.xlabel("Wavelength [nm]", fontweight='bold', fontsize=Font)
# plt.plot(wave2,NormFLUX2*ReNormTerm, 'b.', markersize=1, label='observation UT171108')

########### to get correction term
# plt.plot(wave1,NormFLUX1*(cutResNormFlux[0]/cutResNormFlux1[0]), 'r.', markersize =1)
if np.max(cutResWav2-cutResWav1) != 0:
    sys.exit("We have a problem!!!") 
Correction = (cutResNormFlux2*ReNormTerm)/cutResNormFlux1
Z = np.load('UT170804corrected_spec.npy')
Zz = np.load('corrected_UT170804.npy')
# plt.plot(cutResWav1, cutResNormFlux1*Correction, 'g--')
Gap = np.where((8788.8 < Z[0,:]) & (Z[0,:] < 8892.3))[0]
Z[1,:][Gap] = np.nan
plt.plot(Z[0,:]/10, Z[1,:], 'r.', zorder=4, markersize=1)
WavRng = np.where((wav_rang[0]<wave1) & (wave1<wav_rang[1]))
CorrInterp = np.interp(wave1[WavRng], cutResWav1, Correction)
plt.plot(Zz[0,:]/10, Zz[1,:], 'r--', markersize=1, label='corrected ob1')
# plt.plot(wave1[WavRng], NormFLUX1[WavRng]*CorrInterp, 'r.', markersize=1, label = 'corrected_UT170804')
plt.xticks(fontweight='bold', fontsize=tickFont)
plt.yticks(fontweight='bold', fontsize=tickFont)
plt.xlim((490,940))
plt.legend()

########### to plot correction term
plt.subplot2grid((1,8), (0,3), colspan=2)
CorrInterp 
W= np.load('correct_Term.npy')
WAve, CorrInterp = W[0,:], W[1,:]
plt.plot(WAve/10, CorrInterp)
plt.xlabel("Wavelength [nm]", fontweight='bold', fontsize=Font)
# plt.yticks(ha='left')
# ax.tick_params(axis="y",direction="in", pad=-22)
plt.xticks(fontweight='bold', fontsize=tickFont)
plt.yticks(fontweight='bold', fontsize=tickFont)
plt.xlim((490,940))
plt.subplots_adjust(wspace=0.6, hspace=0)
path = '/Users/chimamcgruder/Research/Papers/WASP96b/ExtractedSpectra/'
plt.savefig(path+'Corr2ndOrderLight.png',bbox_inches='tight')
plt.show()
plt.close()

