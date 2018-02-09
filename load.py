import os
import matlablike as pys
import string
import struct
import re

"""
Basic file loading functionality for Bruker files

"""

"""
EPR Specific file handling
"""
# {{{
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
# }}}

"""
NMR Specific file handling
"""


b0 = r'$B_0$'

#{{{ load the pulse sequence parameters
def load_acqu(filename,whichdim='',return_s = None):
    filename = dirformat(filename)
    if det_type(filename)[0] == 'bruker':
        if return_s is not None:
            return bruker_load_acqu(filename,whichdim=whichdim,return_s = return_s)
        else:
            return bruker_load_acqu(filename,whichdim=whichdim)
    elif det_type(filename)[0] == 'prospa':
        if det_type(filename)[1] == 't1_sub':
            filename = dirformat(filename)
            return prospa_load_acqu(filename+'../')
        else:
            return prospa_load_acqu(filename)
    else:
        raise CustomError(det_type(filename),'is not yet supported')
#}}}

#{{{ Determine type of file Bruker or Prospa 1 or 2 D
## You could turn this into an overall file parser.
def det_type(filename):
    filetype = None
    # WinEPR
    if os.path.exists(filename+'.spc'):
        return ('winepr',True)
    else:
        filename = dirformat(filename)
        files_in_dir = os.listdir(filename)
        # Bruker 2D
        if os.path.exists(filename+'ser'):
            return ('bruker',True)
        # Prospa generic 2D
        elif os.path.exists(filename+'data.2d'):
            return ('prospa',True)
        # specific Prospa formats
        elif any(map((lambda x:'Delay' in x),files_in_dir)):
            return ('prospa','t1')
        elif os.path.exists(filename+'acqu.par'):
            return ('prospa',False)
        elif os.path.exists(filename+'../acqu.par'):
            return ('prospa','t1_sub')
        # Bruker 1D
        elif os.path.exists(filename+'acqus'):
            return ('bruker',False)
        else:
            raise CustomError('WARNING! unidentified file type '+filename)
#}}}

def format_listofexps(args):
    'aux function, just used to decode (filename,listofexpnos) vs. (listoffilenames) in an arbitrary way'
    if type(args[0]) is str: # even if it's just a string, make it a list
        args[0] = [args[0]]
    if len(args) > 1 and ((not isscalar(args[1])) and len(args[1]) == 0): args.pop(1)
    if len(args) > 1:
        if len(args) > 2: raise CustomError('wrong number of args!')
        if isscalar(args[1]): args[1] = [args[1]] # if the second argument is a single file number, make it into a list
        if len(args[0]) > 1: raise CustomError("you can't have both the filename and the expnos be longer than 1")
        filenames = [dirformat(args[0][0]) + '%d'%x for x in args[1]]
    else:
        filenames = args[0]
    return filenames

def load_file(*args,**kwargs):
    'load a file or series of files; load as load_file(filename) or load_file(dirname,expnos)'
    args = list(args)
    #{{{ manually do kwargs
    dimname = '' 
    if 'dimname' in kwargs.keys(): dimname = kwargs['dimname']
    calibration = 1.0
    if 'calibration' in kwargs.keys(): calibration = kwargs['calibration']
    add_sizes = []
    if 'add_sizes' in kwargs.keys(): add_sizes = kwargs['add_sizes']
    add_dims = []
    if 'add_dims' in kwargs.keys(): add_dims = kwargs['add_dims']
    #}}}
    filenames = format_listofexps(args)
    #{{{load all the data into a list
    data = [load_indiv_file(filenames[0],dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    for filename in filenames[1:]:
        data += [load_indiv_file(filename,dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    #}}}
    # for the following, I used to have a condition, but this is incompatible with the pop statement at the end
    newdata = concat(data,dimname) # allocate the size of the indirect array
    #print 'DEBUG concatenated list = ',data
    newdata_shape = ndshape(newdata)
    if all(map((lambda x:det_type(x)[0]=='prospa'),filenames)):
        if hasattr(data[0],'want_to_prospa_decim_correct'):
            if data[0].want_to_prospa_decim_correct is True:
                newdata = prospa_decim_correct(newdata)
    #print 'DEBUG concatenated list before pop = ',data
    if newdata_shape[dimname]==1:
        newdata.popdim(dimname)
    return newdata*calibration

def load_indiv_file(filename,dimname='',return_acq=False,add_sizes=[],add_dims=[]):
    filetype,twod = det_type(filename)
    #{{{ ESR spectra
    if filetype == 'winepr':
        fp = open(filename+'.spc','rb')
        data = fp.read()
        data = array(
                struct.unpack('<%df'%(len(data)/4),data),
                dtype='double')
        v = winepr_load_acqu(filename)
        xpoints = v['RES']
        rg = v['RRG']
        data /= rg
        modulation = v['RMA']
        #data /= modulation
        try:
            data /= v['JNS'] # divide by number of scans
        except:
            pass
        #data /= v['MP'] # divide by power <-- weird, don't do this!
        ypoints = len(data)/xpoints
        if ypoints>1:
            if ypoints != v['REY']:
                raise CustomError('I thought REY was the indirect dim, guess not')
            if dimname=='':
                dimname = v['JEY']
            data = nddata(data,[ypoints,xpoints],[dimname,b0])
        else:
            data = nddata(data,[xpoints],[b0])
        xlabels = linspace(v['HCF']-v['HSW']/2.,v['HCF']+v['HSW']/2.,xpoints)
        if len(data.dimlabels)>1:
            yaxis = r_[0:v['REY']]
            if dimname == 'mw-power-sweep':
                yaxis *= v['MPS']
                yaxis += v['XYLB'] # the starting attenuation
                yaxis = 10**(-yaxis/10.) # convert to linear power
                yaxis *= v['MP']/yaxis[0] # the initial power
                yaxis *= 1e-3 # convert from mW to W
                data.rename('mw-power-sweep','power')
                dimname = 'power'
            data.labels([dimname,b0],[yaxis,xlabels])
            data.reorder([b0,dimname])
        else:
            data.labels([b0],[xlabels])
        return data
    #}}}
    filename = dirformat(filename)
    #{{{ Bruker 2D
    if twod and filetype == 'bruker':
        v = bruker_load_acqu(filename)
        v2 = bruker_load_acqu(filename,whichdim='2')
        td2 = int(v['TD'])
        rg = float(v['RG'])
        td1 = int(v2['TD'])
        td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
        fp = open(filename+'ser','rb')
        data = fp.read()
        if float(v.get('BYTORDA')) == float(1.0):
            ### Read for big endian
            data = array(struct.unpack('>%di'%(len(data)/4),data),
                    dtype='complex128')
        if float(v.get('BYTORDA')) == float(0.0):
            ### Read for little endian
            data = array(struct.unpack('<%di'%(len(data)/4),data),
                    dtype='complex128')
        data = data[0::2]+1j*data[1::2]
        data /= rg
        mydimsizes = [td1,td2_zf/2]
        mydimnames = [dimname]+['t2']
        #print 'DEBUG: data going to nddata =',data
        try:
            data = nddata(data,mydimsizes,mydimnames)
        except:
            size_it_should_be = array(mydimsizes).prod()
            if size_it_should_be > len(data):
                zero_filled_data = zeros(size_it_should_be)
                zero_filled_data[0:len(data)] = data
                data = nddata(zero_filled_data,mydimsizes,mydimnames)
            else:
                new_guess = len(data)/(td2_zf/2)
                print lsafen("WARNING!, chopping the length of the data to fit the specified td1 of ",td1,"points!\n(specified ",zip(mydimnames,mydimsizes),' td2_zf=%d)'%td2_zf)
                #td2_zf_new = 2**ceil(log(td2)/log(2))
                #mydimsizes[1] = td2_zf_new
                #size_it_might_be = array(mydimsizes).prod()
                #print "maybe this works:",size_it_might_be == len(data)
                data = data[0:size_it_should_be]
                data = nddata(data,mydimsizes,mydimnames)
                #raise CustomError("found td1=",td1,"for",filename,"which I don't think is right, because the product of the dimensions",zip(mydimnames,mydimsizes),'=',size_it_should_be,'does not equal the length of the data',len(data),'I think that it should be',len(data)/(td2_zf/2))
        #print 'DEBUG: data straight from nddata =',data
        data = data['t2',0:td2/2] # now, chop out their zero filling
        t2axis = 1./v['SW_h']*r_[1:td2/2+1]
        t1axis = r_[0:td1]
        mylabels = [t1axis]+[t2axis]
        data.labels(mydimnames,mylabels)
        shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the first order phase shift
        data.circshift('t2',shiftpoints)
        data.set_units('t2','s')
        data.set_units('digital')
        data.other_info['title'] = bruker_load_title(filename)
        #print 'DEBUG 2: data from bruker file =',data
        #}}}
        #{{{ Prospa 2D
    elif twod and filetype == 'prospa':
        if twod == 't1_sub':
            v = prospa_load_acqu(filename+'../') # if it's subdirectory format, the file comes from one directory up
            indirect_dim_len = [1]
            indirect_dim_name = [dimname]
            dimshere = 1
        else:
            v = prospa_load_acqu(filename)
            indirect_dim_name = []
            indirect_dim_len = []
            dimshere = 2
        taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3 # this is the t2 dimension, and so is always true
        data = prospa_load_datafile(filename,dims=dimshere)/v['nrScans']#added this 2/20/13 to allow automatic signal averaging
        #{{{ Prospa CPMG
        if v['experiment'].find('cpmg') > -1:
            data = nddata(data,indirect_dim_len+[v['nrEchoes'],v['nrPnts']],indirect_dim_name+['echo','t2'])
            echotime = (r_[0:v['nrEchoes']]+0.5)*v['echoTime']/1e6
            data.labels(indirect_dim_name+['echo','t2'],indirect_dim_len+[echotime,taxis])
            data.want_to_prospa_decim_correct = False
        #}}}
        #{{{ Prospa where 1D subscan is not CPMG
        else:
            data = nddata(data,indirect_dim_len+[v['nrPnts']],indirect_dim_name+['t2'])
            data.labels([dimname,'t2'],[r_[1],taxis])
            data.want_to_prospa_decim_correct = True
        #}}}
        #}}}
        #{{{ bruker 1D
    else:
        if filetype == 'bruker':
            v = bruker_load_acqu(filename)
            td2 = int(v['TD'])
            td1 = 1
            td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
            fp = open(filename+'fid','rb')
            data = fp.read()
            data = array(
                    struct.unpack('>%di'%(len(data)/4),data),
                    dtype='complex128')
            data = data[0::2]+1j*data[1::2]
            rg = v['RG']
            data /= rg
            data = nddata(data,[td1,td2_zf/2],[dimname,'t2'])
            data = data['t2',0:td2/2] # now, chop out their zero filling
            t2axis = 1./v['SW_h']*r_[1:td2/2+1]
            t1axis = r_[1]
            data.labels([dimname,'t2'],[t1axis,t2axis])
            shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the second order phase shift
            #print 'shiftpoints = ',shiftpoints
            data.circshift('t2',shiftpoints)
            # finally, I will probably need to add in the first order phase shift for the decimation --> just translate this
            data.other_info['title'] = bruker_load_title(filename)
        #}}}
        #{{{ prospa 1d
        elif filetype == 'prospa':
            v = prospa_load_acqu(filename)
            data = prospa_load_datafile(filename,dims=1)/v['nrScans']#added this 2/20/13 to allow automatic signal averaging
            data = nddata(data,[v['nrPnts']],['t2'])
            taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
            data.labels(['t2'],[taxis])
        #}}}
        else:
            raise CustomError("can't load this file type $\\rightarrow$ \\verb+%s+"%filename)

    #{{{ return, and if necessary, reorganize
    if len(add_sizes)>0:
        data.labels([dimname],[[]]) # remove the axis, so we can reshape
        #print 'DEBUG: data before chunk = ',data
        data.chunkoff(dimname,add_dims,add_sizes)
        #print 'DEBUG: data after chunk = ',data
        data.labels(add_dims,
                [r_[0:x] for x in add_sizes])
    if return_acq:
        return (data,v,v2)
    else:
        return data
    #}}}

#{{{ t1 axis
def load_t1_axis(file):
    if det_type(file)[0] == 'bruker':
        return bruker_load_t1_axis(file)
    elif det_type(file)[0] == 'prospa':
        return prospa_t1_info(file)[1]
    else:
        raise CustomError('Trying to load T1 axis on a file of unrecognized format!')
#}}}

#{{{ load acqu
def prospa_decim_correct(data):
    #{{{ get rid of the finite rise time    
    data_abs = abs(data)
    otherdims = ndshape(data)
    otherdims.pop('t2')
    for indirect_dim_name in otherdims.dimlabels:
        data_abs.mean(indirect_dim_name)
    data_abs = data_abs.run(argmax,'t2')
    top = int(data_abs.data)
    data.circshift('t2',top)
    #}}}
    print 'Applied prospa decimation correction'
    return data

def prospa_load_acqu(file):
    file = dirformat(file)
    fp = open(file+'acqu.par')
    lines = fp.readlines()
    line_re = re.compile(r'([^ \t]+) *= *(.+)')
    vars = {}
    for j in range(0,len(lines)):
        lines[j] = string.rstrip(lines[j])
        m = line_re.match(lines[j])
        if m:
            exec 'temp = %s'%m.groups()[1]
            vars.update({m.groups()[0]:temp})
        else:
            print "error, acqu.par line not parsed: ",lines[j]
    fp.close()
    return vars
#}}}

#{{{ load_datafile
def prospa_load_datafile(file,dims=1):
    r'''load a prospa datafile into a flat array as a 1D file
    use dims=2 if it's a 2D file'''
    file = dirformat(file)
    if dims == 1:
        fp = open(file+'data.1d','rb')
    elif dims == 2:
        fp = open(file+'data.2d','rb')
    else:
        print 'ERROR: wrong number of dims'
    data = fp.read()
    data = array(struct.unpack('%df'%(len(data)/4),data))
    data = data[7:]
    # the following is junk!!!
    #elif precision=='b':
    #   data = array(struct.unpack('%db'%(len(data)/1),data))
    #   data = data[7*4:]
    #else:
    #   print 'error, precision wrong'
    data = reshape(data,(-1,2))
    data = data[:,0]+1j*data[:,1]
    fp.close()
    return data
#}}}

#{{{ load the data from a t1 file based on file names
def prospa_t1_info(file):
    file = dirformat(file)
    if det_type(file) == ('prospa','t1_sub'):
        file += '../'
    elif not det_type(file) == ('prospa','t1'):
        raise CustomError("You're trying to get prospa T1 info from a file that's not a Prospa T1 file!")
    files = [x for x in os.listdir(file) if os.path.isdir(dirformat(file)+x)]
    file_re = re.compile(r'([0-9]+)msDelay$')
    datafiles = []
    wait_time = []
    print 'DEBUG: prospa is searching for times in the file list',files
    for j in range(0,len(files)):
        m = file_re.match(files[j])
        if m:
            datafiles += [file+files[j]]
            wait_time += [int(m.groups()[0])]
    return datafiles,array(wait_time)*1e-3
#}}}

#{{{ Load T1 axis
def bruker_load_t1_axis(files):
    wait_time = []
    if type(files) is str:
        files = [files]
    #print 'DEBUG: trying to load files: ',files
    for thisfile in files:
        thisfile = dirformat(thisfile)
        thisfiletype = det_type(thisfile)
        if thisfiletype[0] == 'prospa':
            print 'need to copy over code from prospa'
        if thisfiletype[0] == 'bruker':
            wait_time += [bruker_load_vdlist(thisfile)]
        else:
            print 'couldn\'t determine thisfile type'
    return array(wait_time).flatten()
#}}}

### Bruker Specific
#{{{ calculate decimation correction
def bruker_det_phcorr(v):
    if v['DIGMOD']==1:
        gdparray=array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[179,201,533,709,1097,1449,2225,2929,4481,5889,8993,11809,18017,23649,36065,47329,72161,94689,144353,189409,288737],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[184,219,384,602,852,1668,2292,3368,4616,6768,9264,13568,18560,27392,36992,55040,73856,110336,147584,220928,295040]])
        decimarray=array([2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024]) # the -1 is because this is an index, and copied from matlab code!!!
        dspfvs = int(v['DSPFVS'])
        decim = v['DECIM']
        try:
            retval = gdparray[dspfvs,where(decimarray==decim)[0]]/2/decim
        except:
            print r'\begin{tiny}'
            print r'\begin{verbatim}'
            print CustomError('Problem returning',dspfvs,where(decimarray==decim)[0],'from gdparray of size',shape(gdparray),'because decimarray is of size',shape(decimarray))
            print r'\end{verbatim}'
            print r'\end{tiny}'
            retval = 0
        return retval
    else:
        return array([0])
#}}}


def bruker_match_line(line,number_re,string_re,array_re):
    m = number_re.match(line)
    if m:
        retval = (0,m.groups()[0],double(m.groups()[1]))
    else:
        m = string_re.match(line)
        if m:
            retstring = m.groups()[1]
            if retstring[-1]=='>':
                retstring = retstring[:-1]
            retval = (1,m.groups()[0],0,retstring)
        else:
            m = array_re.match(line)
            if m:
                retval = (2,m.groups()[0],(double(m.groups()[1]),double(m.groups()[2])),m.groups()[3])
            else:
                retval = (3,line)
    return retval

#{{{ bruker_load_vdlist
def bruker_load_vdlist(file):
    fp = open(file+'vdlist')
    lines = fp.readlines()
    lines = map(string.rstrip,lines)
    lines = map((lambda x: x.replace('m','e-3')),lines)
    lines = map((lambda x: x.replace('s','')),lines)
    lines = map((lambda x: x.replace('u','e-6')),lines)
    lines = map(double,lines)
    return array(lines)
#}}}

def load_title(file):
    if det_type(file)[0] == 'bruker':
        return bruker_load_title(file)
    else:
        return ''

def bruker_load_title(file):
    file = dirformat(file)
    fp = open(file+'pdata/1/title')
    lines = fp.readlines()
    emptystring = '\r\n'
    while emptystring in lines:
        lines.pop(lines.index(emptystring))
    emptystring = '\n'
    while emptystring in lines:
        lines.pop(lines.index(emptystring))
    return ''.join(lines)

def bruker_load_acqu(file,whichdim='',return_s = True):
    if return_s:
        fp = open(file+'acqu'+whichdim+'s')# this is what I am initially doing, and what works with the matched filtering, etc, as is, but it's actually wrong
    else:
        fp = open(file+'acqu'+whichdim)# this is actually right, but doesn't work with the matched filtering, etc.
    lines = fp.readlines()
    vars = {}
    number_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *([0-9\-\.]+)')
    string_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *<(.*)')
    array_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *\(([0-9]+)\.\.([0-9]+)\)(.*)')
    lines = map(string.rstrip,lines)
    j=0
    retval =  bruker_match_line(lines[j],number_re,string_re,array_re)
    j = j+1
    retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re) #always grab the second line
    while j < len(lines):
        isdata = False
        if retval[0]==1 or retval[0]==2:
            name = retval[1]
            thislen = retval[2]
            data = retval[3]
            while (retval2[0] == 3) and (j<len(lines)): # eat up the following lines
                data += ' '+retval2[1]
                j = j+1
                retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re)
            isdata = True
        elif retval[0]==0:
            name = retval[1]
            data = retval[2]
            isdata = True
        #else:
        #   print 'not a data line:',retval[1]
        if(isdata):
            if retval[0]==2: #if it's an array
                data = data.split(' ')
                if len(data)>0:
                    while '' in data:
                        data.remove('')
                    try:
                        data = map(double,data)
                    except:
                        print "can't map ",data," to float."
                    if len(data)-1!= thislen[1]:
                        print 'error:',len(data)-1,'!=',thislen[1]
            vars.update({name:data})
        # at this point, the string or array data is loaded into data and we have something in retval2 which is definitely a new line
        retval = retval2
        j = j+1
        if j<len(lines):
            retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re)
    fp.close()
    return vars

