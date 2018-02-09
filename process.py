from numpy import *
import matlablike as pys

"""
need to insert pys in appropriate places...

"""


def findPeaksSequential(dataArray, windowSize = 3, threshold = 0.2): # {{{
    """ this finds peaks in the EPR spectrum by a sequential peak finding algorithm. Returns index of peak positions in data set.
    Input:
    dataArray - (array) of EPR spectrum
    windowSize - (int) size of window to find peak in
    threshold - (float) threshold cutoff. Only values above threshold will be recorded as peaks.

    Returns:
    peaks - (list) index values of peak positions
    """

    peaks = []
    counter = 0
    currentMaxIndex = False

    for index in range(len(dataArray)):
        if currentMaxIndex: # We have a currentMaxIndex if we're climbing a peak
            if dataArray[index] >= dataArray[currentMaxIndex]: # the new point is greater than my current maximum point. Reset max index and counter
                currentMaxIndex = index
                counter = 0
            else: # my current max is still the max
                counter += 1
            if counter >= windowSize: # I've found a peak
                peaks.append(currentMaxIndex)
                currentMaxIndex = False
                counter = 0

        else: # We don't have a currentMaxIndex, check to see if we've started climbing again
            if (index > windowSize) and (dataArray[index] >= threshold):
                if dataArray[index] > dataArray[index - 1]:
                    currentMaxIndex = index
    return peaks# }}}

def calcSpinConc(calibrationFile):#{{{
    """
    Use the EPR Double integral value to calculate the sample concentration given a calibration file.
    Format of the calibration file (csv).
    Concentration (uM) X Double Integral Value
    ConcVal              DI Val

    Args:
    CalibrationFile - csv of calibration

    returns:
    calibration - the estimated concentration of the spin system
    """
    openFile = open(calibrationFile,'rt')
    lines = openFile.readlines()
    lines = lines[0].split('\r')
    lines.pop(0)
    concL = []
    diL = []
    for line in lines:
        conc,di = line.split(',')
        concL.append(float(conc))
        diL.append(float(di))
    openFile.close()

    calib = pys.nddata(pys.array(diL)).rename('value','concentration').labels('concentration',pys.array(concL))
    return calib#}}}

def findPeaks(spec,numberOfPeaks,verbose = False):#{{{
    """
    DEPRECIATED! 

    Find the position of the peaks and valleys of the EPR spectrum given the number of peaks to look for. 
    The function returns the total peak to peak width of the spectrum, given more than one peak, as well as the center field and linewidth.

    args:
    spec - an nddata set of the EPR spectrum. The EPR spectrum should be the data and the field values should be placed in an axis named 'field'
    numberOfPeaks - an integer. The number of peaks to find, for nitroxide this should be 3.

    """
    print "'findPeaks()' DEPRECIATED use 'findPeaksSequential()'"
    peaks = []
    valleys = []
    hrf = linspace(spec.getaxis('field').min(),spec.getaxis('field').max(),10000)
    smash = spec.copy().interp('field',hrf).runcopy(real) # use an interpolated higher res spec to get a more accurate esitmate of linewidth
    hrs = smash.copy()
    #smash -= average(spec.data)
    for i in range(numberOfPeaks): 
        peak = smash.data.argmax()
        peaks.append(peak)
        valley = smash.data.argmin()
        valleys.append(valley)
        # remove from peak
        #find the high bound
        notCrossed=True
        count = 0
        dimSize = len(smash.data)
        while notCrossed:
            if peak + count <= 0:
                lowBound = peak+count
                notCrossed = False
            else:
                if float(smash['field',peak+count].data) <= 0.0:
                    lowBound = peak+count
                    notCrossed = False
            count-=1
        # find the low bound
        notCrossed=True
        count=0
        while notCrossed:
            if peak + count >= dimSize: # check to make sure you haven't wandered off the spectrum
                highBound = peak+count
                notCrossed = False
            else:
                if float(smash['field',peak+count].data) <= 0.0:
                    highBound = peak+count
                    notCrossed = False
            count+=1
        smash['field',lowBound:highBound] = 0.0

        # remove from valley
        #find the high bound
        notCrossed=True
        count = 0
        while notCrossed:
            if valley + count <= 0:
                lowBound = valley+count
                notCrossed = False
            else:
                if float(smash['field',valley+count].data) >= 0.0:
                    lowBound = valley+count
                    notCrossed = False
            count-=1
        # find the low bound
        notCrossed=True
        count=0
        while notCrossed:
            if valley + count >= dimSize: # check to make sure you haven't wandered off the spectrum
                highBound = valley+count
                notCrossed = False
            else:
                if float(smash['field',valley+count].data) >= 0.0:
                    highBound = valley+count
                    notCrossed = False
            count+=1
        smash['field',lowBound:highBound] = 0.0
        if verbose:
            pys.plot(smash)
    peak = pys.nddata(hrs.data[peaks]).rename('value','field').labels('field',hrs.getaxis('field')[peaks])
    valley = pys.nddata(hrs.data[valleys]).rename('value','field').labels('field',hrs.getaxis('field')[valleys])
    # Calculate relevant parameters
    peak.sort('field')
    valley.sort('field')
    return peak,valley
#}}}

def workupCwEpr(eprName,spectralWidthMultiplier = 1.25,numPeaks=3,EPRCalFile=False,firstFigure=[]): #{{{ EPR Workup stuff
    """
    This is a beast of a function...
    Perform the epr baseline correction and double integration.

    Args:
    eprName - string - full name of the EPR file without the file extension.
    spectralWidthMultiplier - a multiplier to determine the integration width based on the width of the spectrum (difference between lowField and highField peak)
    numPeaks - the number of peaks of the spectrum
    EPRCalFile - a calibration file in csv format
    firstFigure - a figure list instance (see fornotebook)

    Returns:
    spec - nddata - the EPR spectra with other info set to the EPR params dict.
    lineWidths - list - the EPR linewidths
    spectralWidth - double - the EPR peak to peak spectral width
    centerField - double - the centerfield
    doubleIntZC - nddata - the double integral spectrum
    doubleIntC3 - nddata - the double integral spectrum with cubic baseline correction
    diValue     - float  - the double integral value
    spinConc    - float  - the spin concentration from the double integral. Must have an EPRCalFile
    amplitudes  - array  - amplitude of each peak in the spectrum
    """
    firstFigure.append({'print_string':r'\subparagraph{EPR Spectra %s}'%eprName + '\n\n'})
    eprFileName = eprName.split('\\')[-1]
    # Pull the specs, Find peaks, valleys, and calculate things with the EPR spectrum.#{{{
    spec,normalized = returnEPRSpec(eprName,resample=True)
    peak,valley = findPeaks(spec,numPeaks)
    lineWidths = valley.getaxis('field') - peak.getaxis('field') 
    amplitudes = peak.data - valley.data
    spectralWidth = peak.getaxis('field').max() - peak.getaxis('field').min() 
    # determine the center field
    if numPeaks == 2:
        centerField = (peak.getaxis('field')[0] + lineWidths[0]/2. + peak.getaxis('field')[1] + lineWidths[1]/2.)/2.
    elif numPeaks == 3:
        centerField = peak.getaxis('field')[1] + lineWidths[1]/2.
    specStart = centerField - spectralWidthMultiplier*spectralWidth
    specStop = centerField + spectralWidthMultiplier*spectralWidth
    print "\nI calculate the spectral width to be: ",spectralWidth," G \n"
    print "I calculate the center field to be: ",centerField," G \n"
    print "I set spectral bounds of: ", specStart," and ", specStop," G \n"#}}}
    if normalized == 'bad':
        print "The spectra is not normalized by the receiver gain"

    # Baseline correct the spectrum #{{{
    baseline1 = spec['field',lambda x: x < specStart].copy()
    baseline2 = spec['field',lambda x: x > specStop].copy()
    specBase = array(list(baseline1.data) + list(baseline2.data))
    fieldBase = array(list(baseline1.getaxis('field')) + list(baseline2.getaxis('field')))
    specBase = pys.nddata(specBase).rename('value','field').labels('field',fieldBase)
    ### Calculate 0th, 1st, and 3rd order baselines
    baseline = average(specBase.data)
    spec.data -= baseline # zeroth order correct the spectrum

    # Plot the results
    firstFigure = pys.nextfigure(firstFigure,'EPRSpectra')
    pys.plot(spec,'m',alpha=0.6)
    pys.plot(peak,'ro',markersize=10)
    pys.plot(valley,'ro',markersize=10)
    pys.plot(spec['field',lambda x: logical_and(x>specStart,x<specStop)],'b')
    pys.title('Integration Window')
    pys.ylabel('Spectral Intensity')
    pys.xlabel('Field (G)')
    pys.giveSpace(spaceVal=0.001)
    #}}}

    ### Take the first integral #{{{
    absorption = spec.copy().integrate('field')#}}}

    # Fit the bounds of the absorption spec to a line and subtract from absorption spectrum.#{{{
    baseline1 = absorption['field',lambda x: x < specStart]
    baseline2 = absorption['field',lambda x: x > specStop]
    fieldBaseline = array(list(baseline1.getaxis('field')) + list(baseline2.getaxis('field')))
    baseline = pys.concat([baseline1,baseline2],'field')
    baseline.labels('field',fieldBaseline)
    # Do the first order correction
    c1,fit1 = baseline.polyfit('field',order = 1)
    fit1 = pys.nddata(array(c1[0] + absorption.getaxis('field')*c1[1])).rename('value','field').labels('field',absorption.getaxis('field'))
    correctedAbs1st = absorption.runcopy(real) - fit1.runcopy(real)
    c3,fit3 = baseline.polyfit('field',order = 3)
    fit3 = pys.nddata(array(c3[0] + absorption.getaxis('field')*c3[1] + (absorption.getaxis('field')**2)*c3[2] + (absorption.getaxis('field')**3)*c3[3])).rename('value','field').labels('field',absorption.getaxis('field'))
    correctedAbs3rd = absorption.runcopy(real) - fit3.runcopy(real)
    #}}}

    # Set the values of absorption spec outside of int window to zero.#{{{
    zeroCorr = correctedAbs1st.copy()
    zeroCorr['field',lambda x: x < specStart] = 0.0
    zeroCorr['field',lambda x: x > specStop] = 0.0#}}}

    # Plot absorption results#{{{
    firstFigure = pys.nextfigure(firstFigure,'Absorption')
    pys.plot(absorption,label='uncorrected')
    pys.plot(fit1,label='1st order fit')
    pys.plot(fit3,label='3rd order fit')
    pys.plot(correctedAbs1st,label='1st corr')
    pys.plot(correctedAbs3rd,label='3rd corr')
    pys.plot(zeroCorr,label='zero cut')
    pys.title('Absorption Spectrum')
    pys.ylabel('Absorptive Signal')
    pys.xlabel('Field (G)')
    pys.giveSpace(spaceVal=0.001)
    pys.legend()
    #}}}

    # Calculate and plot the double integral for the various corrections you've made #{{{
    doubleInt = absorption.copy().integrate('field')
    doubleIntC1 = correctedAbs1st.copy().integrate('field')
    doubleIntC3 = correctedAbs3rd.copy().integrate('field')
    doubleIntZC = zeroCorr.copy().integrate('field')
    diValue = doubleIntC3.data.max()
    print "\nI calculate the double integral to be: %0.2f\n"%diValue

    firstFigure = pys.nextfigure(firstFigure,'DoubleIntegral')
    pys.plot(doubleInt,label='uncorrected')
    pys.plot(doubleIntC1,label='1st corrected')
    pys.plot(doubleIntC3,label='3rd corrected')
    pys.plot(doubleIntZC,label='zero corrected')
    pys.legend(loc=2)
    pys.title('Double Integral Results')
    pys.ylabel('Second Integral (arb)')
    pys.xlabel('Field (G)')
    pys.giveSpace(spaceVal=0.001)
    #}}}
    
    # If the calibration file is present use that to calculate spin concentration#{{{
    if normalized == 'good':
        if EPRCalFile:
            calib = calcSpinConc(EPRCalFile)
            ### Fit the series and calculate concentration
            c,fit = calib.polyfit('concentration')
            spinConc = (diValue - c[0])/c[1]
            # Plotting 
            firstFigure = pys.nextfigure(firstFigure,'SpinConcentration')
            pys.plot(calib,'r.',markersize = 15)
            pys.plot(fit,'g')
            pys.plot(spinConc,diValue,'b.',markersize=20)
            pys.title('Estimated Spin Concentration')
            pys.xlabel('Spin Concentration')
            pys.ylabel('Double Integral')
            ax = pys.gca()
            ax.text(spinConc,diValue - (0.2*diValue),'%0.2f uM'%spinConc,color='blue',fontsize=15)
            pys.giveSpace()
        else:
            spinConc = None
            #}}}
    return spec,lineWidths,spectralWidth,centerField,doubleIntZC,doubleIntC3,diValue,spinConc,amplitudes
    #}}}

#{{{ optimize the zero-order phase the right way
def phaseopt(curve):
    curve = curve.copy()
    #{{{ find bestindex once
    phases = linspace(-pi/2,pi/2,100).reshape(1,-1) # this should work w/out altering the sign
    rotated_data = (curve.reshape(-1,1))*exp(-1j*phases)
    success = (real(rotated_data)**2).sum(axis=0)/((imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
    bestindex = where(success==max(success))[0][0]
    #}}}
    #{{{ find the bestindex based on that 
    if (bestindex>phases.shape[1]-2):
        bestindex = phases.shape[1]-2
    phases = linspace(
            phases[0,bestindex-1],
            phases[0,bestindex+1],
            100).reshape(1,-1)
    rotated_data = (curve.reshape(-1,1))*exp(-1j*phases)
    success = (real(rotated_data)**2).sum(axis=0)/((imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
    bestindex = where(success==max(success))[0][0]
    #}}}
    return exp(-1j*phases[0,bestindex])
#}}}

#{{{ integrate --> new integration function
def integrate(file,expno,
        integration_width=1e3,
        dimname = 'power',
        intpoints=None,
        show_image=True,
        filter_edge = 10,
        center_peak=True,
        use_baseline=True,
        plot_check_baseline=False,
        filter_direct = False,
        return_noise=True,
        show_integral = False,
        indiv_phase = False,
        scale_plot = False,
        peak_within = 1e3,
        bandpass = None,
        abs_image = False,
        max_drift = 1e3,
        first_figure = None,
        pdfstring = '',
        phnum = [],
        phchannel = [],
        returnIntData= False,
        offset_corr = 0,
        forceGlitch=False,
        timeZeroGlitch = True,
        breakIntegrate = False,
        test_drift_limit = False):
    r'''new integration function, which replaces integrate_emax, and is used to integrate data, as for Emax and T1 curves'''
    #print lsafen("DEBUG: yes, integrate was called")
    figurelist = figlistini_old(first_figure)
    if type(plot_check_baseline) is bool:
        if plot_check_baseline:
            plot_check_baseline = 0 
        else:
            plot_check_baseline = -1
    phcycdims = ['phcyc%d'%j for j in range(1,len(phnum)+1)]
    data = load_file(file,expno,dimname = dimname, add_sizes = phnum,add_dims = phcycdims) # load the data
    #{{{ offset correction
    if breakIntegrate:
        return data,figurelist
    if type(offset_corr) is list:
        offset_corr = array(offset_corr)
    if type(offset_corr) is ndarray:
        data.data -= data['t2',offset_corr].copy().mean('t2').mean(dimname).data
    elif offset_corr > 0: # number of points to use for digitizer offset (zero glitch) correction
        offset_verbose = False
        if offset_verbose:
            nextfigure(figurelist,'beforecorr')
            image(abs(data))
            colorbar()
            print 'DEBUG: attempting offset correction by',offset_corr
        if array(offset_corr).dtype is dtype('float64'):
            data.data -= data['t2',lambda x: x > offset_corr].copy().mean('t2').mean(dimname).data
            if offset_verbose: print 'which is a double'
        else:
            data.data -= data['t2',-offset_corr:].copy().mean('t2').mean(dimname).data
        if offset_verbose:
            print '\n\n'
            nextfigure(figurelist,'aftercorr')
            image(abs(data))
            colorbar()
    #}}}
    # see_fid obsolete by rg_check
    # also remove all debug statements
    #print 'DEBUG: before phcyc, figlist is',lsafen(figurelist)
    data,figurelist = phcyc(data,names = phcycdims,selections = phchannel, show_plot = ['t2',(lambda x:abs(x)<peak_within)],first_figure = figurelist,pdfstring = pdfstring,bandpass = bandpass) # ft along t2, applying phase cycle where necessary
    if timeZeroGlitch:
        # Remove the zero glitch that occurs before the actual FID.
        data.ift('t2',shift = True)
        ### Pull out the initial receiver glitch and zero fill then end of the data set.
        zeroglitch = data.runcopy(abs)
        zeroglitch.sum(dimname)
        if forceGlitch:
            glitch = forceGlitch
        else:
            glitch = zeroglitch.run(argmax,'t2').data
            print "I remove %0.1f points in the beginning of the spectrum"%glitch
        dataList = []
        for count in range(len(data.getaxis(dimname))):
            zeroFree = data[dimname,count]
            zeroFree = nddata(array(list(zeroFree['t2',glitch:].data) + [complex(0.0)]*glitch)).rename('value','t2').labels('t2',data.getaxis('t2'))
            dataList.append(zeroFree)
        data = concat(dataList,dimname).labels(dimname,data.getaxis(dimname))
        data.reorder('t2',dimname)
        figurelist = nextfigure(figurelist,'zeroGlitchRemoval' + pdfstring)
        plot(data)
        title('After cutting and zero filling')
        data.ft('t2',shift = True)

    data_shape = ndshape(data) # this is used to shape the output
    #{{{ abs w/ ma SNR, so we can pick the peak
    data_abs = data.copy()
    data_abs['t2',(lambda x: abs(x)<peak_within)].data *= 0
    #{{{ apply the matched filter to maximize our SNR while picking the peak
    data_abs.ift('t2',shiftornot=True) # ft along t2
    filter = matched_filter(data_abs,'t2',decay_rate = 3)
    data_abs *= filter
    data_abs.ft('t2',shiftornot=True) # ft along t2
    #}}}
    data_abs = abs(data_abs)
    #}}}
    #{{{ generate topavg --> the index at the top of the average
    data_mean = data_abs.copy()
    data_mean.mean(dimname) # note that we are taking the mean over the abs here, which would not be great for the noise, but the center should still be in the right place
    data_mean.run(argmax,'t2') # put the index of the top peak there
    topavg = int32(data_mean.data)
    #}}}
    #{{{ generate center --> an array with the value at the center of each scan
    data_center = data_abs.copy() # since data_abs is currently not used, but I want to use it to do matched filtered integration, really need to make a separate variable here
    f = data_center.getaxis('t2')
    # If the spectrum is off resonance you need to account for this somehow...
    # you could try setting max drift to the average resonant position, but you need to drop any outliers.
    data_center['t2',abs(f-f[topavg])>max_drift] = 0 # we need to keep the indeces in place, but don't want to pick anything too far out of the way
    if test_drift_limit:
        figurelist = nextfigure(figurelist,'driftLimitTest' + pdfstring)
        plot(data_center.reorder(['t2',dimname]))
        xlim(f[topavg]-max_drift,f[topavg]+max_drift)
        title('You Should See Every Peak In this Window')
    data_center_sum = data_center.copy()
    data_center_sum.sum_nopop('t2')
    data_center.data[:] /= data_center_sum.data # now the it sums to one so we can get a weighted average over the indeces
    f_ind = r_[0.0:double(size(f))].reshape(data_center.getaxisshape('t2')) # make a list of the indeces along the frequency axis
    #data_center.data[:] *= f_ind # multiply by that list
    #data_center.sum('t2') # sum, so that we would return the mean of the list if the spectrum were flat
    data_center.run(argmax,'t2')
    #center = int32(round(data_center.data))
    center = int32(array(map(round,data_center.data)))
    #}}}
    #{{{ if integration points unspec, pull out the integration width
    if intpoints==None:
        df = data.getaxis('t2').copy()
        df = df[1]-df[0]
        intpoints = floor(integration_width/(df))
    #}}}
    data_shape['t2'] = intpoints*2+1
    newdata = []
    newnoise = []
    center[center<intpoints] = intpoints # prevent a bug where the integration range exceeds the spectrum
    # plotting here removed, since this is done by phcyc
    #{{{baseline correction
    if use_baseline:
        #data.data += 0.1 #for debug baseline
        for j in range(0,ndshape(data)[dimname]):
            if plot_check_baseline == j: # if I've passed this index, show this baseline
                baseline_data = baseline_spectrum(data[dimname,j].copy(),center[j],intpoints,showplots = True) # call baseline_spectrum with center and intpoints, which should already be defined
                error_plot('Wanted to check baseline on scan ',j)
            else:
                baseline_data = baseline_spectrum(data[dimname,j].copy(),center[j],intpoints) # call baseline_spectrum with center and intpoints, which should already be defined
                data[dimname,j].data[:] -= baseline_data.data.flatten()
                if any(isnan(data[dimname,j].data[:])):
                    print 'isnan!!'
                if any(isinf(data[dimname,j].data[:])):
                    print 'isinf!!'
    #}}}
    #{{{ actually pull the integral points and the standard deviation of the noise
    #print 'DEBUG: center points are',center
    plot_noise = [] # array where we can put the noise for plotting
    other_info_retval = data.other_info # hack to compensate for the fact that other info isn't carried through the slice
    for j in range(0,data_shape[dimname]):
        newdata += [data[dimname,j,'t2',center[j]-intpoints:center[j]+intpoints+1]]
        if return_noise:
            #newnoise += [data[dimname,j,'t2',10:10+intpoints]]
            #{{{ grab intpoints to the left of the spectrum 
            temp = center[j]+r_[0:intpoints]-2*intpoints
            temp = temp[temp>0]
            #}}}
            #{{{ grab intpoints to the right of the spectrum 
            list_of_noise_indeces = center[j]+1+r_[0:intpoints]+intpoints
            list_of_noise_indeces = list_of_noise_indeces[list_of_noise_indeces<ndshape(data)['t2']]
            #}}}
            list_of_noise_indeces = int32(r_[temp,list_of_noise_indeces])
            #{{{ pull the noise data and calculate the standard deviation of the noise
            #print 'DEBUG: shape of data',ndshape(data),j,list_of_noise_indeces
            temp = data[dimname,j,'t2',list_of_noise_indeces]
            if show_image:
                plot_noise += [temp.copy()]
            temp.data = abs(temp.data-temp.data.mean()) # need to explicitly do the abs, since the data is complex
            temp.data **= 2
            temp.mean('t2') # I DO NOT take the sqrt here, because it's taken at the very end
            #}}}
            newnoise += [temp.copy()] # find the standard deviation of the noise which we have pulled --> it should be independent of the number of points that we're using
    newdata = concat(newdata,dimname)
    newdata.other_info = other_info_retval # hack to compensate for the fact that other info isn't carried through the slice
    if return_noise:
        newnoise = concat(newnoise,dimname)
        if show_image:
            plot_noise = concat(plot_noise,dimname)
    #}}}
    if show_image:
        plot_newdata = newdata.copy() # make a backup for plotting
    #{{{ autophase
    if not indiv_phase:
        phaseoptval = phaseopt(newdata.data)
        newdata.data *= phaseoptval
        if show_image:
            plot_newdata.data *= phaseoptval
            if return_noise:
                plot_noise.data *= phaseoptval
        # do NOT rotate the noise data to be returned, since it's real!
    else:
        for j in range(0,len(newdata.data)):
            phcorr =  newdata[dimname,j]
            phcorr /= abs(phcorr)
            try:
                if show_image:
                    plot_newdata[dimname,j] *= phcorr
                newdata[dimname,j] *= phcorr
                if return_noise:
                    #newnoise[dimname,j] *= phcorr # again, don't rotate the real noise
                    if show_image:
                        plot_noise[dimname,j] *= phcorr
            except:
                print 'shape of newdatacopy',ndshape(newdatacopy)
                print 'shape of newdata',ndshape(newdata)
    newdata.sum('t2') # integrate --> note that I converted this to a sum!
    #}}}
    #{{{ show what we're integrating
    if show_image:
        figurelist = nextfigure(figurelist,'intpeaks' + pdfstring)
        #print "DEBUG intpeaks figurelist =",figurelist,"gcf = ",gcf().number
        plot_newdata.reorder(['t2',dimname])
        def maybescale(x):
            if scale_plot:
                return x/newdata
            else:
                return x
        if scale_plot:
            plot_newdata /= newdata
            plot_noise /= newdata
        plot(plot_newdata,alpha=0.5)
        title('Peaks, zoomed in to integration region')
        if return_noise:
            plot_noise.reorder(['t2',dimname])
            plot_color_counter(0)
            plot(plot_noise['t2',0:intpoints],'-',alpha=0.1)
            plot_color_counter(0)
            plot(plot_noise['t2',intpoints:],'-',alpha=0.1)
        if show_integral:
            #{{{this does work to plot the integral
            plot_newdata.integrate('t2') # this is apparently a function to do integral with all the correct bells and whistles
            #gridandtick(gca())
            ax = gca()
            myxlim = copy(ax.get_xlim())
            twinx()
            plot(plot_newdata,'k',alpha=0.1)
            gca().set_xlim(myxlim)
            #}}}
    #}}}
    #print lsafen("DEBUG: ready to return from integrate")
    #{{{ return integral and the noise of the 
    if return_noise:
        number_of_integral_points = 2.*intpoints+1. # the integral is a sum over this many points
        newdata.set_error(sqrt(newnoise.data.flatten()*number_of_integral_points))
    if first_figure == None:
        return newdata
    elif returnIntData:
        return newdata,figurelist,plot_newdata
    else:
        return newdata,figurelist
    #}}}
#}}}

#{{{ deal with manual phase cycling
def phcyc(data,names=[],selections=[],remove_zeroglitch=None,show_plot = False,first_figure = None,pdfstring = '',bandpass = None):
    figurelist = figlistini_old(first_figure)
    data.ft('t2',shift=True)
    if (bandpass is not None):
        data = data['t2',lambda x: abs(x)<bandpass]
    if len(names)>0:
        data.ft(names)
        if remove_zeroglitch == None:
            remove_zeroglitch = True
    if remove_zeroglitch:
        index = argmin(abs(data.getaxis('t2')-0)) # need to incorporate this functionality into the abstracted function indexer
        indexlist = ['t2',index]
        for j,name in enumerate(names):
            indexlist += [name,0]
        data[tuple(indexlist)] = 0
    if show_plot:
        nextfigure(figurelist,'phcycchan' + pdfstring)
        allindexes = list(data.dimlabels)
        for name in names:
            if name in allindexes:
                allindexes.remove(name)
                allindexes = [name]+allindexes
        allindexes.remove('t2')
        data.reorder(allindexes+['t2'])
        if show_plot == True: # if passed something other than just "true", use that to subscript data
            image(data)
        else:
            image(data[tuple(show_plot)])
        titlestr = 'Raw data'
        if len(names)>0:
            titlestr += ' by phase cycle channel'
        title(titlestr)
    other_info_retval = data.other_info
    for j,name in enumerate(names):
        data = data[name,selections[j]]
    data.other_info = other_info_retval # this and the line that saves other_info_retval above are a hack because nddata doesn't do this correctly
    if first_figure == None:
        return data
    else:
        return data,figurelist
#}}}

#{{{ process_t1
def process_t1(file,expno,usebaseline = None,showimage = None,plotcheckbaseline = None,saturation = False,first_figure = None,pdfstring = '',t1_offset_corr = None,verbose = False,showlinear = False,**kwargs):
    """
    DEPRECIATED ?
    No you could utilize this in the ODNP workup software
    """
    #{{{ hack it, since it only actually takes a single file 
    file = format_listofexps([file,expno])
    if len(file) > 1: raise CustomError('I don\'t think this can handle more than one file at a time')
    expno = []
    #}}}
    figurelist = figlistini_old(first_figure)
    #{{{ legacy kwargs
    if showimage != None:
        show_image = showimage
    if usebaseline != None:
        use_baseline = usebaseline
    #}}}
    if type(file) is str:
        file = [file]
    titlestr = load_title(file[0])
    wait_time = load_t1_axis(file[0])
    if t1_offset_corr is not None:
        kwargs.update({'offset_corr':t1_offset_corr})
    integral,figurelist = integrate(file,expno,first_figure = figurelist,pdfstring = pdfstring,**kwargs)
    t1name = r'$t_1$'
    integral.rename('power',t1name)
    integral = t1curve(integral,fit_axis = t1name) # make this into an integral class, which fits along the dimension t1
    if ndshape(integral)[t1name] < len(wait_time):
        print '\n\nNote: ',t1name,'axis shorter than list of delays'
        wait_time = wait_time[0:ndshape(integral)[t1name]]
    integral.labels([t1name],[wait_time]) # before, I had to sort them manually, but now, I don't
    if verbose: print 'DEBUG wait times:',integral.getaxis(t1name)
    integral.sort(t1name)
    #{{{ finally, show the fit  
    figurelist = nextfigure(figurelist,'t1'+pdfstring)
    taxis = wait_time
    integral.data *= phaseopt(integral.data)
    plot(integral.runcopy(real),'ko')
    plot(integral.runcopy(imag),'yo')
    integral.makereal() # otherwise, it won't fit
    integral.fit()
    if verbose: print 'DEBUG: after fit, fit coeff is',integral.fit_coeff
    plot(integral.eval(300)) # evaluate the fit function on the axis taxis
    #{{{ for now, do not plot the modified versions
    #plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*1.2],taxis),'y')
    #plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*0.8],taxis),'y')
    #}}}
    ax = gca()
    text(0.5,0.75,integral.latex(),transform = ax.transAxes,size = 'x-large', horizontalalignment = 'center',color = 'r')
    title(titlestr)
    #}}}
    if showlinear:
        #{{{ and the straight line plot
        figurelist = nextfigure(figurelist,'t1straight'+pdfstring)
        #print 'linear data:',integral.linear().data
        plot(integral.linear(),'o')
        plot(integral.linear(taxis))
        #print ndshape(integral.linear())
        #}}}
    if first_figure == None:
        return integral # there is never a return_fit, since the fit is stored in the class itsself
    else:
        return integral,figurelist # there is never a return_fit, since the fit is stored in the class itsself
#}}}

#{{{ process cpmg 
def process_cpmg(file,dimname=''):
    """
    DEPRECIATED ?
    No I don't think so. Just not currently used
    r'this just takes an FT along T2 and pulls the signal at the frequency where it\'s maxed'
    """
    data = load_file(file,dimname=dimname)
    data.ft('t2')
    findmax = abs(data)
    findmax.mean_all_but(['t2'])
    findmax = findmax.run(argmax,'t2').data
    data = data['t2',findmax]
    data.mean_all_but(['echo','t1'])
    data.data *= phaseopt(data.data) # I added this in, not sure why it was gone!
    return data
#}}}

#{{{ matched filter
def matched_filter(data,along_dim,decay_rate = 1,return_fit=False):
    r'''take ift'd data, and apply the matched filter to it
    This takes time domain data, fits it to decaying exponent and applys the fit to the data to apodize the signal
    
    '''
    #{{{ actually find the filter   
    data_abs = abs(data)
    timeaxis = data_abs.getaxis('t2')
    labels = list(data.dimlabels)
    labels.pop(labels.index(along_dim))
    for thisdim in labels:
        data_abs.mean(thisdim)
    p = exp_fit(timeaxis,data_abs.data)
    #}}}
    #{{{ actually apply the filter, note this does not actually apply the filter..
    filter = ndshape(data)
    for thisdim in labels:
        filter.pop(thisdim)
    filter = filter.alloc()
    if (not return_fit): # don't mess with it if we want to check the fit
        p[1] /= decay_rate
    filter.data = exp_fitfunc(p,data.getaxis(along_dim).copy())
    #print 'for matched filter, the x axis is ',data.getaxis(along_dim).copy()
    if not return_fit:
        filter.data[:] -= filter.data.flatten()[-1] # drop so that end is at zero (since we have a noise baseline)
        filter.data[:] /= filter.data.flatten()[0] # normalize
    filter.labels(['t2'],[data_abs.getaxis('t2')])
    return filter
    #}}}
#}}}

def __baseline_gen_L(data):# {{{
    x = data.getaxis('t2').copy()
    x_norm = max(abs(x))
    x /= x_norm # normalize, otherwise we get ridiculously large higher order terms
    L = array([ones(shape(x)),x,x**2/2,x**3/6,x**4/24,x**5/120]).T
    #L = array([ones(shape(x)),x,x**2/2]).T
    return x,x_norm,L# }}}

def baseline_spectrum_peakpick(data,showplots=False,threshold=10,check_filter=False,set_errorlevel=False):# {{{
    #print 'diagnose: start baseline'
    data = data.copy()
    prelim_offset = mean(r_[data.data[0],data.data[1],data.data[-2],data.data[-1]])
    data.data[:] -= prelim_offset
    data.ift('t2',shiftornot=True)
    if check_filter:
        clf()
        #plot(abs(data))
        plot(matched_filter(data,'t2',return_fit=True))
        legend(['fits to']) #legend(['data','fits to'])
        twinx()
        plot(matched_filter(data,'t2'))
        legend(['filter with'])
        return
    #{{{ make the abs of a broadened spectrum   
    data_widefilt = data.copy() * matched_filter(data,'t2',decay_rate=10)
    data_widefilt.ft('t2',shiftornot=True)
    #data_widefilt.data *= phaseopt(data_widefilt.data)
    #data_widefilt.data = real(data_widefilt.data)
    #data_widefilt.data /= sign(data_widefilt.data[argmax(abs(data_widefilt.data))])
    data_widefilt.data -= (data_widefilt.data[0]+data_widefilt.data[-1])/2
    data_widefilt = abs(data_widefilt)
    mask = (data_widefilt.data<(data_widefilt.data.max()/threshold)) # generate the mask according to the threshold
    #}}}
    #{{{ ft the data
    data.ft('t2',shiftornot=True)
    #}}}
    if sum(mask)==0:
        erroronnopeak = False
        if erroronnopeak:
            legendstring = []
            plot(abs(data))
            legendstring += ['mf data']
            plot(data_widefilt)
            legendstring += ['wide filter']
            legend(legendstring)
            error_plot("Error -- not able to find any non-baseline data")
        else:
            print "Warning, fit entire spectrum to baseline.\n\n"
        mask = (data_widefilt.data<(data_widefilt.data.max()/1.5))
        mask[:] = True
    data_baseline = data.copy()
    data_baseline.data = data_baseline.data[mask]
    data_baseline.axis_coords[data_baseline.dimlabels.index('t2')] = data_baseline.getaxis('t2')[mask]
    legendstring = []
    x,x_norm,L = __baseline_gen_L(data_baseline) # return a normalized x axis, the value used to normalize it, and the array of normalized polynomials
    try:
        fit_coeff = dot(pinv(L,rcond=1e-5),data_baseline.data) # L * fit_coeff = data
    except:
        raise CustomError(maprep('problem inverting:',shape(L),shape(data_baseline.data)))
    #print 'diagnose: inverted'
    if set_errorlevel:
        if any(abs(dot(L,fit_coeff))>set_errorlevel):
            showplots = True
    if showplots:
        #plot(abs(data))
        #legendstring += ['data']
        plot(data_widefilt)
        legendstring += ['wide filter']
        plot(abs(data_baseline))
        legendstring += ['baseline portion']
        show_L = False
        if show_L:
            plot(x*x_norm,L)
            legendstring += map(
                    (lambda x:'L'+str(x)),
                    range(1,1+L.shape[1])
                    )
        plot(x*x_norm,abs(dot(L,fit_coeff)))
        legendstring += ['fit to baseline']
    x,x_norm,L = __baseline_gen_L(data)
    if showplots:
        plot(x*x_norm,dot(L,fit_coeff))
        legendstring += ['entire fit']
        legend(legendstring,'best')
    baseline_data = nddata(prelim_offset+dot(L,fit_coeff),[size(x),1],['t2','power'])
    baseline_data.labels(['t2'],[x])
    #print 'diagnose: shape of baseline ',ndshape(baseline_data)
    return baseline_data# }}}

def baseline_spectrum(data,center,points,showplots=False):# {{{
    data = data.copy()
    #{{{ here, I should define a mask in the same way I do in integrate, just inverted, so that I DON'T use the spectrum
    mask = bool8(zeros(shape(data.data)))
    mask[center-points:center+points+1] = True # this is an exact copy of the indeces used to create newdata in integrate
    mask = ~mask
    #}}}
    #{{{ splice out data_baseline --> the data which makes up the baseline
    data_baseline = data.copy()
    data_baseline.data = data_baseline.data[mask]
    data_baseline.axis_coords[data_baseline.dimlabels.index('t2')] = data_baseline.getaxis('t2')[mask]
    #}}}
    #{{{ perform the leastsquare fit
    x,x_norm,L = __baseline_gen_L(data_baseline) # return a normalized x axis, the value used to normalize it, and the array of normalized polynomials
    try:
        fit_coeff_Re = dot(pinv(L,rcond=1e-5),real(data_baseline.data)) # L * fit_coeff_Re = real(data)
    except:
        raise CustomError(maprep('problem inverting:',shape(L),shape(data_baseline.data)))
    try:
        fit_coeff_Im = dot(pinv(L,rcond=1e-5),imag(data_baseline.data)) # L * fit_coeff_Im = imag(data)
    except:
        raise CustomError(maprep('problem inverting:',shape(L),shape(data_baseline.data)))
    fit_coeff = fit_coeff_Re + 1j * fit_coeff_Im
    #}}}
    # here, i deleted set_errorlevel, since I can't remember what it does, so it must not be important
    #{{{ generate the matrices that span the full dataset and show the plots with all the info
    if showplots:
        clf() # this is in case I'm not running in the notebook, and want a decent error plot
        legendstring = []
        plot(abs(data),alpha=0.5)
        legendstring += ['data']
        plot(abs(data_baseline))
        legendstring += ['baseline portion']
        show_L = False
        if show_L:
            plot(x*x_norm,L)
            legendstring += map(
                    (lambda x:'L'+str(x)),
                    range(1,1+L.shape[1])
                    )
        plot(x * x_norm,abs(dot(L,fit_coeff)))
        legendstring += ['fit to baseline']
    x,x_norm,L = __baseline_gen_L(data)
    if showplots:
        plot(x*x_norm,dot(L,fit_coeff))
        legendstring += ['entire fit']
        data_forplot = data.copy()
    #}}}
    #{{{ generate and return the baseline curve
    baseline_data = nddata(dot(L,fit_coeff),[size(x),1],['t2','power'])
    baseline_data.labels(['t2'],[x])
    #}}}
    #{{{ show what the resulting data should look like
    if showplots:
        data.data[:] -= baseline_data.data.flatten()
        plot(abs(data))
        legendstring += ['baseline corrected data']
        legend(legendstring,'best')
    #}}}
    return baseline_data# }}}

#{{{ plot_noise
def plot_noise(path,j,calibration,mask_start,mask_stop,rgmin=0,k_B = None,smoothing = False, both = False, T = 293.0,plottype = 'semilogy',retplot = False):
    '''plot noise scan as resistance'''
    filename = r'%s%d'%(path,j)
    try:
        data = load_file(filename,calibration=calibration)
    except:
        raise CustomError('error loading file'+filename)
    k_B = 1.3806504e-23
    data.ft('t2',shift = True)
    newt2 = r'F2 / $Hz$'
    data.rename('t2',newt2)
    v = bruker_load_acqu(r'%s%d/'%(path,j))
    dw = 1/v['SW_h']
    dwov = dw/v['DECIM']
    rg = v['RG']
    de = v['DE']
    aq = v['TD']*dw
    if rg>rgmin:
        plotdata = abs(data)
        plotdata.data **= 2
        johnson_factor = 4.0*k_B*T
        plotdata.data /= (aq*johnson_factor)
        t = data.getaxis(newt2)
        mask = logical_and(t>mask_start,
            t<mask_stop)
        try:
            avg = plotdata.data[mask].mean() 
        except IndexError:
            raise CustomError('error trying to mask for the average because the mask is',mask,'of shape',shape(mask),'for shape(plotdata)=',shape(plotdata.data))
        retval = []
        if both or not smoothing:
            pval = plot(plotdata,'-',alpha=0.5,plottype = plottype)
            retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+'$t_{dw}$ %0.1f $t_{dwov}$ %0.1f RG %d, DE %0.2f, mean %0.1f'%(dw*1e6,dwov*1e6,rg,de,avg)]
            axis('tight')
        if smoothing:
            # begin convolution
            originalt = plotdata.getaxis(newt2).copy()
            plotdata.ft(newt2,shift = True)
            sigma = smoothing
            siginv = 0.5*sigma**2 # here, sigma is given in the original units (i.e. what we're convolving)
            t = plotdata.getaxis(newt2)
            g = exp(-siginv*t.copy()**2) # we use unnormalized kernel (1 at 0), which is not what I thought!
            plotdata.data *= g
            plotdata.ift(newt2,shift = True)
            t = plotdata.getaxis(newt2).copy()
            t[:] = originalt
            # end convolution
            pval = plot(plotdata,'-',alpha=0.5,plottype = plottype)
            retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, DE %0.2f, mean %0.1f'%(dwov*1e6,rg,de,avg)]
            axis('tight')
        if retplot:
            return pval,retval
        else:
            return retval
    else:
        return []
#}}}

