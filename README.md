# 2nd_Order_Light_Correction
Code used to model and correct for second-order light that contaminants stellar spectroscopic images on IMACS

A spectrograph differentiates all the wavelengths of white light via diffraction. That diffraction occurs multiple times, each time with reduced intensity, called 'orders'. To observe exoplanet atmospheres, we observe minute changes in the host stars' stellar spectra due to the exoplanets’ atmosphere. Specifically, we used the Inamori Magellan Areal Camera and Spectrograph (IMACS) mounted on the Baade Magellan telescope to observe the stellar spectra. 

For our observations, we only use the first-order diffracted light and initially thought the higher-order light was diffracted off our detectors. We soon discovered that this was not the case and we needed to use a light-blocking filter to prevent second-order light contamination. However, to use data that was obtained before we included the blocking filter, I had to model out the light contamination. To do this, I created a script that modeled the light contamination, this script.

This script:
1) First uses iSpec to model the normalized continuum of the spectra observed. This is needed so the general spectroscopic light curve can be obtained without spectroscopic features found in specific observations. So things like changing spectroscopic lines don’t change the correction.
2) Then the normalized continuum of the older unfiltered spectra is convolved with the throughput profile of the filter that was used on newer observations. This produces the spectral profile of the observation if the filter had been added, but with any second-order light still present. 
3) Given the only difference between the newer observations *with* the filter and the older observation convolved with the filter (step 2) is the second-order light, I used the ratio of the two continua as a correction continuum for the second-order light. 
4) This correction can then be applied to all data observed without the blocking filter


Scripts:
Corr2ndOrderLight - the main script to correct second-order light contamination \n
AnalyzeDat.py - the majority of this script is used for other analyses. However, ‘degrade_resolution()’ function is used here \n
2nd_Order_light - a plot comparing a test observation with the filter (UT171108), without the filter (UT170804), and without the filter but convolved with the throughput of the filter. Comparing the spectrum of UT171108 (blue) and UT170804 after being convolved with the filter’s throughput (green) highlights the effect of the second order light. \n

Dependencies:
astropy - https://www.astropy.org/ \n
iSpec - https://www.blancocuaresma.com/s/iSpec \n
