import os
import matlablike as pys

"""
Basic file loading functionality for Bruker files

"""

def loadEPRFits(fileName):# {{{
    """
    Load in the fitting data that the multi component NLSL fitting program returns.

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

def loadEPRExpDict(fileName,extension = '.par',verbose=False):#{{{
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

def loadEPRExpDictDSC(fileName):#{{{
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

def loadFTEseemTrace(fileName):#{{{
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

def loadPulseEPRSpec(fileName,expType = 'fieldDomain',spect='Elexsys'):#{{{
    """
    This actually opens any single dimension pulse experiment performed on the Elexsys or Specman and returns an nddata.

    You're writing this because bruker's file format changes too much. You also included support for specman 1-D filetype

    spect - string - the spectrometer name i.e. 'Elexsys' or 'Specman'

    """
    if spect == 'Elexsys':
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
    elif spect == 'Specman':
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

def loadCwEPRSpec(fileName,doNormalize = True, resample=False): #{{{
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

def loadEPRt2TwoDim(path,name,extension='.DTA',runsToCut=False,firstFigure=[],showPlots=True):# {{{
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

