"""
This script will calculate the double integral of a derivative EPR spectrum that is normalized by both the number of scans and the receiver gain which gives reproducible values.

You should make this dump the spectrum and the double integral to the database with a searchable sample number. This way you could calculate all ODNP + EPR information by looking for you sample.

"""

import matlablike as pys
from numpy import *
import csv
import os

pys.close('all')

# Various Definitions and classes#{{{
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

def loadEPRFits(fileName):# {{{
    """
    Load in the fitting data that the multi component fitting program returns.

    """
    fileHandle = open(fileName,'r')
    lines = fileHandle.readlines()
# find out how many data lists I need.
    numSets = len(lines[0].split('\r\n')[0].split('\t'))
# the structure is C0 - field, C1 - data, C2 - fit result, C3 - weights, C4 - component 1, C5 - component 2, C6 and on are more components. 
    numComps = numSets - 4
    toStore = zeros((len(lines), numSets))
    for count, line in enumerate(lines):
        line = line.split('\r\n')[0].split('\t')
        for count1, item in enumerate(line):
            toStore[count, count1] = item
    rawData = pys.nddata(toStore[:,1]).rename('value','field').labels('field',toStore[:,0])
    fit = pys.nddata(toStore[:,2]).rename('value','field').labels('field',toStore[:,0])
    components = {}
    for i in range(numComps):
        components.update({'%i'%i: pys.nddata(toStore[:,i+4]).rename('value','field').labels('field',toStore[:,0])})
    return rawData, fit, components# }}}

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

def returnEPRExpDict(fileName,extension = '.par',verbose=False):#{{{
    """
    Return all of the experiment parameters stored in the '.par' file output by the Bruker

    Args:
    fileName - string - full file name from top directory.

    Returns:
    expDict - dictionary - Keys are keys from bruker par files, values are everything else matched to the corresponding key.
    """
    openFile = open(fileName + extension,'r') # read the par
    lines = openFile.readlines()
    expDict = {}
    if extension == '.par':
        lines = lines[0].split('\r')
    for line in lines:
        try:
            if verbose:
                print "Debug: ",line
            if '=' in line:
                print "pulling for eq"
                splitData = line.split('=')
                key = splitData[0].split(' ')[0]
                value = splitData[1]
            else:
                splitData = line.split(' ')
                key = splitData.pop(0)
                value = splitData.pop(0)
                for data in splitData:
                    value += data
            expDict.update({key:value})
        except:
            pass
    return expDict#}}}

def returnEPRExpDictDSC(fileName):#{{{
    """
    This returns the exp dict stored in the dsc files written by xepr
    """
    openFile = open(fileName + '.DSC','r') # read the par
    lines = openFile.readlines()
    expDict = {}
    for count,line in enumerate(lines):
        cut = line.split('\n')[0]
        if 'end defs' in cut:
            break
        if '\t' in cut:
            try:
                key,value = cut.split('\t')
                expDict.update({key.strip():value.strip()})
            except:
                pass
        elif '=' in cut:
            cut =  cut.split('=')
            key = cut[0].strip()
            value = cut[1].split(';')[0].strip()
            expDict.update({key:value})
        else:
            splits = cut.split(' ')
            key = splits.pop(0).strip()
            value = splits
            expDict.update({key:value})
    return expDict#}}}

def returnFTEseemTrace(fileName):#{{{
    """
    Opens the FT ESEEM trace that is saved in the PEPP Matlab processing program.

    returns an nddata with frequency dimension in MHz
    """
    openFile = open(fileName,'r+')
    lines = openFile.readlines()
    freq = []
    signal = []
    for line in lines:
        line = filter(None,line.split('\n')[0].split(' '))
        freq.append(float(line[0]))
        signal.append(float(line[1]))
    signal = pys.nddata(array(signal)).rename('value','MHz').labels('MHz',array(freq))
    signal = signal['MHz',lambda x: x >= 0.0]
    return signal#}}}

def returnPulseEPRSpec(fileName,expType = 'fieldDomain',spect='xepr'):#{{{
    """
    This actually opens any single dimension pulse experiment performed on the Elexsys and returns an nddata.

    You're writing this because bruker's file format changes too much. You also included support for specman 1-D filetype

    spect - string - the spectrometer name i.e. 'xepr' or 'specMan'

    """
    if spect == 'xepr':
        expDict = returnEPRExpDictDSC(fileName)
        specData = fromfile(fileName+'.DTA','>d') # or if it is a DTA file read that instead
        real = []
        imaginary = []
        for i in range(0,len(specData),2):
            real.append(specData[i])
            imaginary.append(specData[i+1])

        data = array(real)+1j*array(imaginary)
        numScans = int(expDict.get('n'))
        numAcc = int(expDict.get('h'))
        if expType == 'fieldDomain':
            centerField = float(expDict.get('CenterField')[-2])
            sweepWidth = float(expDict.get('SweepWidth')[-2])
            spec = pys.nddata(data).rename('value','field').labels('field',linspace(centerField-sweepWidth/2,centerField+sweepWidth/2,len(data)))
        elif expType == 'timeDomain':
            timeStart = float(expDict.get('XMIN'))
            timeStop = float(expDict.get('XWID'))
            unit = str(expDict.get('XUNI'))
            if unit == "'ns'":
                multiplier = 1e-9
            elif unit == "'us'":
                multiplier = 1e-6
            else:
                multiplier = 1
            spec = pys.nddata(data).rename('value','time').labels('time',linspace(timeStart,timeStop,len(data)))
            spec.other_info.update({'timeUnit':unit})

        spec /= numScans
        spec /= numAcc
        return spec
    elif spect == 'specMan':
        openFile = open(os.path.abspath(fileName),'r+')
        lines = openFile.readlines()
        lines.pop(0)
        time = []
        real = []
        imag = []
        for line in lines:
            line = filter(None,line.split('\n')[0].split(' '))
            time.append(float(line[0]))
            real.append(float(line[1]))
            imag.append(float(line[2]))
        spec = pys.nddata(array(real)+1j*array(imag)).rename('value','time').labels('time',array(time))
        return spec

#}}}

def returnEPRSpec(fileName,doNormalize = True, resample=False): #{{{
    """ 
    *** This code is crappy

    Right now you try to incorporate stuff for xepr cw scans and you do it in a try except loop which is not the way to do this!!! This is bad code. Fix when you aren't in a rush!
    # you might want to force a choice on spc or dta so that you can tell the necessary workup to perform for the given file type as the normalization is different.

    ***

    Return the cw-EPR derivative spectrum from the spc and par files output by the winEPR program.
    If doNormalize is set to True (recommended) this will normalize the spectral values to the number of scans run as well as the receiver gain settings. This is a more reproducible value as is independent of the two settings which may vary.

    args:

    fileName - (sting) full file name not including extension. e.g. '/Users/StupidRobot/exp_data/ryan_cnsi/epr/150525_ConcentrationSeries/200uM_4OHT_14-7mm'

    returns: 

    1-D nddata dimensioned by field values of spectrum, and containing the EPR experimental parameters as other_info.
    """
    # Open the spc and par files and pull the data and relevant parameters
    #try:
    expDict = returnEPRExpDict(fileName)
    specData = fromfile(fileName+'.spc','<f') # read the spc
    sizeY = expDict.get('SSY')
    xU = 'field'
    if sizeY: # this is a two dimensional data set
        sizeY = int(sizeY)
        sizeX = int(expDict.get('SSX'))
        yU = expDict.get('XYUN')
        specData = specData.reshape((sizeY,sizeX))
    if expDict.get('HCF'):
        centerSet = float(expDict.get('HCF'))
    elif expDict.get('XXLB'):
        lowBound = float(expDict.get('XXLB'))
        width = float(expDict.get('XXWI'))
        centerSet = lowBound + width/2.
    else:
        centerSet = float(expDict.get('GST'))
            
    sweepWidth = float(expDict.get('HSW'))
    if doNormalize:
        numScans = expDict.get('JNS') # I'm not sure if this is right
        if numScans:
            numScans = float(numScans)
        else:
            numScans = 1
        specData /= numScans # normalize by number of scans
        if expDict.get('RRG'):
            rg = float(expDict.get('RRG'))
            modAmp = float(expDict.get('RMA'))
            specData /= modAmp # normalize by modulation amplitude
            specData /= rg # normalize by receiver gain
            normalized = 'good'
        else:
            normalized = 'bad'
    else:
        normalized = 'None'
    #except:
    #    expDict = returnEPRExpDictDSC(fileName)
    #    specData = fromfile(fileName+'.DTA','>c') # or if it is a DTA file read that instead
    #    centerSet = float(expDict.get('CenterField').split(' ')[0])
    #    sweepWidth = float(expDict.get('SweepWidth').split(' ')[0])
    #    numScans = float(expDict.get('NbScansAcc')) # Yea bruker just changes things...
    #    rg = float(expDict.get('RCAG'))
    #    if doNormalize:
    #        specData /= rg
    #    normalized = 'good'
    #    sizeY = False

    # calculate the field values and normalize by the number of scans and the receiver gain and return an nddata
    # The data is two dimensional so find second dimension and 
    if sizeY:
        fieldVals = pys.r_[centerSet-sweepWidth/2.:centerSet+sweepWidth/2.:sizeX*1j]
        LB = float(expDict.get('XYLB'))
        width = float(expDict.get('XYWI'))
        yDim = pys.r_[LB : LB + width : sizeY*1j]
        if yU == 'dB': # Change it to power mW.
            yDim = 197.9 * 10**(-1*yDim / 10)
            yU = 'mW'

        dataShape = pys.ndshape([sizeY, sizeX],[yU, xU])
        data = dataShape.alloc(dtype='float')
        data.data = specData
        spec = data
        spec.labels([yU, xU],[yDim, fieldVals])
    else:
        fieldVals = pys.r_[centerSet-sweepWidth/2.:centerSet+sweepWidth/2.:len(specData)*1j]
        spec = pys.nddata(specData).rename('value',xU).labels(xU,fieldVals)
    if resample:
        # down sample the data to 512. This is for output to the multicomponent fitting program.
        newField = pys.r_[centerSet-sweepWidth/2.:centerSet+sweepWidth/2.:512*1j]
        spec = spec.interp(xU,newField)
    spec.other_info = expDict
    return spec,normalized #}}}

def findPeaks(spec,numberOfPeaks,verbose = False):#{{{
    """
    Find the position of the peaks and valleys of the EPR spectrum given the number of peaks to look for. 
    The function returns the total peak to peak width of the spectrum, given more than one peak, as well as the center field and linewidth.

    args:
    spec - an nddata set of the EPR spectrum. The EPR spectrum should be the data and the field values should be placed in an axis named 'field'
    numberOfPeaks - an integer. The number of peaks to find, for nitroxide this should be 3.

    """
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

def returnt2TwoDim(path,name,extension='.DTA',runsToCut=False,firstFigure=[],showPlots=True):# {{{
    fileName = path+name
    # check for the time file
    if os.path.isfile(fileName+'Time.npy'):
        time = load(fileName+'Time.npy')*1e-9
    else: # pull data from the DSC file
        expDict = returnEPRExpDictDSC(fileName)
        start = float(expDict.get('d1'))*1e-9
        step = float(expDict.get('d30'))*1e-9
        xLen = int(expDict.get('XPTS'))
        yLen = int(expDict.get('YPTS'))
        time = r_[start:start + step * xLen: xLen * 1j]

    if extension == '.DTA':
        # grab data and dump everything into an nddata
        dataShape = pys.ndshape([xLen,yLen],['time','run'])
        data2d = dataShape.alloc(dtype='complex')
        data2d.labels(['time','run'],[time,r_[0:yLen]])
        specData = fromfile(fileName+extension,'>d') # or if it is a DTA file read that instead
        dataList = []
        for count in arange(0,len(specData),2):
            dataList.append(specData[count]+1j*specData[count+1])
        for dim in range(yLen):
            data = array(dataList[dim * xLen: (dim + 1) * xLen])
            data2d['run',dim] = data
    elif extension == '.npy':
        specData = load(fileName+extension)
        yLen,xLen=shape(specData)
        dataShape = pys.ndshape([xLen,yLen],['time','run'])
        data2d = dataShape.alloc(dtype='complex')
        data2d.labels(['time','run'],[time,r_[0:yLen]])
        for dim in range(yLen):
            data2d['run',dim] = specData[dim]
        
    if showPlots:
        firstFigure = pys.nextfigure(firstFigure,'AccumEchoDecayCurvesMag' + name)
        pys.image(data2d.runcopy(abs))
        pys.title('Magnitude Relaxation')
        R,T = pys.meshgrid(data2d.getaxis('run'),data2d.getaxis('time'))
        firstFigure = pys.nextfigure(firstFigure,'AccumEchoDecayCurves' + name)
        CS = pys.contour(T,data2d.data,R,len(data2d.getaxis('run')),alpha = 0.2)
        pys.xlabel('time')
        pys.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='on') # labels along the bottom edge are off
        pys.ylabel('$Magnetization$')
        pys.title('Total Scans')
        pys.colorbar()
        if runsToCut:
            data2d = data2d['run',lambda x: x > runsToCut]
            firstFigure = pys.nextfigure(firstFigure,'AccumEchoDecayCurvesCut' + name)
            R,T = pys.meshgrid(data2d.getaxis('run'),data2d.getaxis('time'))
            CS = pys.contour(T,data2d.data,R,len(data2d.getaxis('run')),alpha = 0.2)
            pys.xlabel('time')
            pys.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='on',      # ticks along the bottom edge are off
                top='off',         # ticks along the top edge are off
                labelbottom='on') # labels along the bottom edge are off
            pys.ylabel('$Magnetization$')
            pys.colorbar()
            pys.title('Run Selection')
    return data2d# }}}
#}}}

### Import the files - for now this is hard coded and this only works with ASCII files, you need to change this so you can use the par files as well.
eprPath = '/Users/StupidRobot/exp_data/ryan_cnsi/epr/150707_CheY_8MUreaSeries/'
eprName = 'M17C_8MUrea_10-9mm'
eprName = eprPath + eprName

def workupCwEpr(eprName,spectralWidthMultiplier = 1.25,numPeaks=3,EPRCalFile=False,firstFigure=[]): #{{{ EPR Workup stuff
    """
    Perform the epr baseline correction and double integration.

    Args:
    eprName - string - full name of the EPR file without the file extension.

    Returns:
    spec - nddata - the EPR spectra with other info set to the EPR params dict.
    lineWidths - list - the EPR linewidths
    spectralWidth - double - the EPR peak to peak spectral width
    centerField - double - the centerfield
    doubleIntZC - nddata - the double integral spectrum
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

#spec,lineWidths,spectralWidth,centerField,doubleIntZC,doubleIntC3,diValue,spinConc,amplitudes = workupCwEpr(eprName)


