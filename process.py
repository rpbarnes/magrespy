from numpy import *
import matlablike as pys


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
