# just make this a library of all NMR reading software
from matlablike import *
from nmrfit import *
import re
import string
import struct
import os
import fornotebook
from scipy.io import loadmat
import datetime
import time

def OUTPUT_notebook():
    return True
#{{{ general, non file-format specific functions
#{{{ auto_steps
def show_powerlog(filename,threshold = -35,extra_t1_problem = False,figure_list = None,**kwargs):
    """
    DEPRECIATED 
    """
    #{{{ generate the powers for the T1 series
    figure_list = fornotebook.figlistini(figure_list)
    #figure_list.text('\n\nCheck the $T_1$ powers:\n\n')
    figure_list.next('powerlog_raw')
    figure_list.next('checkpowerlog')
    figure_list.next('powerlog')
    t1_dbm,figure_list.figurelist = auto_steps(filename,
        threshold = threshold,
        first_figure = figure_list.figurelist,# note that auto_steps is old-style, so I pull the list out of the figlist class
        **kwargs)
    t1mask = ones(len(t1_dbm),dtype='bool')
    if extra_t1_problem == True:
        t1mask[-1] = 0 # this is the line you sometimes want to uncomment
    figure_list = check_autosteps(threshold,t1_dbm,figure_list = figure_list,mask = t1mask)
    figure_list.text('\n\nresulting array:'+lsafen(t1_dbm))
    #}}}
    return figure_list,t1_dbm,t1mask

def check_autosteps(threshold,values,figure_list = None,mask = None):
    """
    DEPRECIATED 
    """
    if figure_list == None:
        figure_list = figlistl()
    figure_list.next('checkpowerlog')
    x = r_[r_[0:len(values)],r_[0:len(values)]+0.99]
    y = r_[values,values]
    s = argsort(x)
    plot(x[s],y[s],'b',linewidth=1)
    if mask is not None:
        if sum(mask) > 0:
            x = r_[r_[0:len(values)][~mask],r_[0:len(values)][~mask]+0.99,r_[0:len(values)][~mask]+0.991]
            y = r_[values[~mask],values[~mask],NaN*ones(shape(values[~mask]))]
            s = argsort(x)
            plot(x[s],y[s],'r',linewidth=2)
    expand_y()
    expand_x()
    ax = gca()
    ylims = list(ax.get_ylim())
    ylims[0] = threshold
    ax.set_ylim(ylims)
    gridandtick(ax)
    ylims = list(ax.get_ylim())
    ylims[0] = threshold
    ax.set_ylim(ylims)
    title('Interpretation:\n(red values are not used)')
    xlabel('experiment number')
    ylabel('power (dBm)')
    return figure_list
def auto_steps(filename,threshold = -35, upper_threshold = 30.0, t_minlength = 0.5*60,minstdev = 0.1,showplots = True, showdebug = False,t_start=0,t_stop=60*1000,tolerance = 2,t_maxlen = inf,return_lastspike = False,first_figure = None):
    r'Plot the raw power output in figure 1, and chop into different powers in plot 2'
    figurelist = figlistini_old(first_figure)
    v = loadmat(filename)
    p_ini = v['powerlist']
    t_ini = v['timelist']
    #{{{ plot the raw power log, then just pull out the times we're interested in
    if showplots:
        figurelist = nextfigure(figurelist,'powerlog_raw')
        plot(t_ini/60,p_ini)
        ax = gca()
        ylim = list(ax.get_ylim())
        ylim[0] = -35
        ax.set_ylim(ylim)
        gridandtick(ax)
        title(filename)
    mask = t_ini < t_stop
    minsteps = sum(t_ini<t_minlength)
    maxsteps = sum(t_ini<t_maxlen)
    mask = logical_and(mask,p_ini < upper_threshold)
    mask = logical_and(mask,t_ini > t_start)
    mask = logical_and(mask,isfinite(p_ini))
    mask = logical_and(mask,~(p_ini==-40.0))
    p_ini = p_ini[mask]
    t_ini = t_ini[mask]
    #}}}
    #a = whereblocks(logical_or(p_ini>threshold,~(isfinite(p_ini))))
    a = whereblocks(logical_or(p_ini>threshold,~(isfinite(p_ini))))
    a = a[argmax(map(len,a))]
    t = t_ini[a]
    p = p_ini[a]
    flattened = NaN*zeros(len(p)+minsteps)
    flattenedstd = NaN*zeros(len(p)+minsteps)
    track_curmean = NaN*zeros(len(p)+minsteps)
    track_threshold = NaN*zeros(len(p)+minsteps)
    maxp = max(p)
    minp = min(p)
    spikes = minp*ones(len(p)+minsteps)
    plotdev = zeros(len(p)+minsteps)
    plotstd = zeros(len(p)+minsteps)
    powerlist = []
    stdlist = []
    nextpos = 0
    while nextpos < len(t)-2:
        nextpos += 1 # just to skip a couple points so I don't grab the rise 
        figurelist = nextfigure(figurelist,'powerlog')
        # grab the stdev and average for the minimum number of steps
        blockstart = nextpos
        subset= p[nextpos:minsteps+nextpos]
        curmean = mean(subset[isfinite(subset)])
        stdev = std(subset[isfinite(subset)])
        if stdev < minstdev:
            stdev = minstdev
        # now that we have a decent stdev, just test to see that every step is less than the stdev
        track_curmean[nextpos] = curmean
        nextpos += minsteps
        #{{{ iterate over blocks
        while (nextpos < len(t)-1) and (abs(p[nextpos]-curmean)<tolerance*stdev or ~isfinite(p[nextpos])):
            subset= p[blockstart:nextpos]
            curmean = mean(subset[isfinite(subset)])
            stdev = std(subset[isfinite(subset)])
            if stdev < minstdev:
                stdev = minstdev
            if showdebug:
                plotdev[nextpos] = p[nextpos]-curmean
                plotstd[nextpos] = tolerance*stdev 
            track_curmean[nextpos] = curmean
            track_threshold[nextpos] = tolerance*stdev
            nextpos += 1
            ##{{{ if we get up to the maximum length, I need to break twice the maximum length into two jumps
            if t[nextpos]-t[blockstart] > t_maxlen:
                if t_minlength > t_maxlen:
                    raise CustomError('for auto_steps, minlength can\'t be greater than maxlen')
                #print 'DEBUG: maxlen triggered, nextpos=',nextpos
                biggestjump = blockstart+minsteps+1+argmax(abs(diff(t[blockstart+minsteps:blockstart+2*maxsteps]))) # find the biggest jump over the interval of two lengths
                if biggestjump-blockstart < 2*minsteps: # if I can't make two steps before the biggest jump
                    nextpos = blockstart + minsteps + 1 + argmax(abs(diff(t[blockstart+minsteps:blockstart+maxsteps])))
                    break
                if t[biggestjump]-t[blockstart] > t_maxlen: # otherwise, the biggest jump happens after another jump
                    #print 'DEBUG: greater than'
                    try:
                        nextpos = blockstart + minsteps + 1 + argmax(abs(diff(t[blockstart+minsteps:biggestjump-minsteps])))
                    except:
                        #figlisterr(figurelist,basename = pdfstring)
                        raise CustomError("I don't have room to make two minimum length steps of ",minsteps,"between the start of the block at ",blockstart," and the biggest jump at",biggestjump)
                    break
                else: # and sometimes the biggest jump happens before another jump, but longer than twice the minlen
                    #print 'DEBUG: less than'
                    nextpos = biggestjump
                    break
            ##}}}
        #}}}
        #{{{ need to recalculate the mean
        subset= p[blockstart:nextpos]
        curmean = mean(subset[isfinite(subset)])
        #}}}
        try:
            track_curmean[nextpos] = curmean
        except:
            raise CustomError('track\_curmean is only',len(track_curmean),'and t is only',len(t),'but nextpos is',nextpos)
        #uptohere = flattened[0:nextpos]
        #uptohere[uptohere==0] = cursum/curnum
        #flattened[nextpos-1] = cursum/curnum
        subset = flattened[0:nextpos] # subset of all the flattened stuff up to here
        subset[isnan(subset)] = curmean
        spikes[nextpos] = maxp # show a spike to see clearly where the data is divided
        try:
            lastspike = t[nextpos-minsteps]-t[-1]
        except:
            raise CustomError('len(spikes)=',len(spikes),'len(t)+minsteps=',len(t)+minsteps)
        subset = flattenedstd[0:nextpos] # subset of all the flattened stuff up to here
        subset[isnan(subset)] = stdev
        powerlist += [curmean]
        stdlist += [stdev]
        #flattened[nextpos-1] = curmean
    if showplots:
        plot(t/60,spikes[0:len(t)],'r')
        plot(t_ini/60,p_ini.flatten(),'x-b')
        plot(t/60,flattened[0:len(t)]+2*flattenedstd[0:len(t)],'g',alpha=0.5)
        plot(t/60,track_curmean[0:len(t)],'y',alpha=0.5)
        plot(t/60,track_curmean[0:len(t)]-track_threshold[0:len(t)],'k',alpha=0.5)
        plot(t/60,track_curmean[0:len(t)]+track_threshold[0:len(t)],'k',alpha=0.5)
        plot(t/60,flattened[0:len(t)]-2*flattenedstd[0:len(t)],'g',alpha=0.5)
        title(filename)
    if showdebug:
        if showplots:
            twinx()
        plot(t/60,plotdev,'k')
        plot(t/60,plotstd,'r')
        title('Power meter log')
    gridandtick(gca())
    retval = [array(powerlist)]
    if return_lastspike == True:
        retval += [lastspike]
    if first_figure is not None:
        retval += [figurelist]
    if len(retval) > 1:
        return tuple(retval)
    else:
        return retval[0]
#}}}

#}}} This is modified to just take power and time arrays instead of a power file, Depending on the shiz this may not be necessary, power_nddata should be an nddata with an axis labeled 't'


#}}}

#}}}
#{{{ wrappers/generic functions to load acq and data files
#{{{ load an nddata structure for a 2d set -- give the data needed to load
#}}}
#}}}
#{{{ lower level functions
#{{{ routines specific to prospa
#}}}
#{{{ routines specific to Bruker
#{{{ load acqu

#{{{ bruker_load_acqu
#}}}
#{{{ winepr_load_acqu
def winepr_load_acqu(file):
    fp = open(file+'.par','rU') # the U automatically converts dos format
    lines = fp.readlines()
    vars = {}
    line_re = re.compile(r'([_A-Za-z0-9]+) +(.*)')
    lines = map(string.rstrip,lines)
    #lines = [j.rstrip('\n') for j in lines] # because it's just \n, even on windows
    v = {'DRS':4096,'RES':1024,'HSW':50}
    for line in lines:
        m = line_re.match(line)
        if m is None:
            raise RuntimeError('Warning:',lsafen(repr(line)),'does not appear to be a valid WinEPR format line, and I suspect this is a problem with the terminators!')
        else:
            name = m.groups()[0]
            value = m.groups()[1]
            try:
                value = int(value)
            except:
                try:
                    value = double(value)
                except:
                    pass
            v[name]=value
    jss = long(v['JSS'])
    parameters = [ 'DUAL', '2D', 'FT', 'MAN0', 'MAN1', 'PROT', 'VEPR', 'POW', 'ABS', 'FTX', 'FTY', 'POW2', 'ABS2']
    parameters = map((lambda x: 's_'+x),parameters)
    masks = [ 0x00000001L, 0x00000002L, 0x00000004L, 0x00000008L, 0x00000010L, 0x00000020L, 0x00000040L, 0x00000080L, 0x00000100L, 0x00000200L, 0x00000400L, 0x00000800L, 0x00001000L]
    values = map((lambda x: x&jss),masks)
    values = map(bool,values)
    values = map(bool,values)
    v.update(dict(zip(parameters,values)))
    return v
#}}}
#}}}
#}}}
#}}}
#{{{ higher level functions
# standard_epr moved to h5nmr
#{{{ fitting functions

#}}}
#{{{ routines for processing emax and t1 type curves
#{{{ rg_check --> check the receiver gain
def rg_check(file,expno,number_of_samples = 75,dynamic_range = 19,first_figure = None,show_complex = True):
    """
    DEPRECIATED ?
    r'''show the fid, plotted vs. the max possible value to check the max of the receiver gain in figure 1
    in figure 2, plot in the complex plain to check digitization limit'''
    """
    fl = figlistini_old(first_figure)
    # the following was load_emax, and I have not yet checked it
    data = load_file(file,expno,dimname = 'power',printinfo = False) # load the data
    listoffiles = format_listofexps([file,expno])
    params = load_acqu(listoffiles[0]) # load the parameters
    maxfloat = 3.4028234e38
    maxint = 2147483648
    max_for_this_dynamic_range= 2**(dynamic_range-1) # the -1 is the bit used for the sign
    data.reorder(['t2','power'])
    spec_type = det_type(listoffiles[0])[0]
    if spec_type == 'bruker':
        data.data *= params['RG'] # scale back up by RG, so we get raw numbers
    data.data /= max_for_this_dynamic_range # scale down by the max int number, so we test digitization
    nextfigure(fl,'rg1')
    if normalize == True:
        data['t2',:] /= data['t1',1:2]
    plot(data['t2',0:number_of_samples],'k')
    plot(data['t2',0:number_of_samples]*(-1j),'b',alpha=0.3)
    title('receiver gain upper value check')
    axis('tight')
    if show_complex:
        nextfigure(fl,'rg2')
        OLDplot(real(data['t2',0:number_of_samples].data),imag(data['t2',0:number_of_samples].data),'.',markersize=0.6)
        #polar(angle(data['t2',0:number_of_samples].data),abs(data['t2',0:number_of_samples].data),'.',markersize=0.1)
        title('receiver gain minimum value check')
        axis('tight')
    if first_figure == None:
        return
    else:
        return fl
#}}}
#}}}
#{{{ load the data from a emax series based on input array
def print_info(filename,also = {}):
    """
    DEPRECIATED ?
    'print the info for a file: to add other parameters, call with the "also" dictionary, with the description as a key, and the variable name or variable name, index number pair as the key'
    """
    filetype,twod = det_type(filename)
    if filetype == 'bruker':
        v = bruker_load_acqu(dirformat(filename))
        f = open(dirformat(filename)+'pulseprogram','r')
        ppginfo = f.read()
        f.close()
        for m in re.finditer(r'\b([pd])([0-9]+)\b',ppginfo):
            also.update({m.group():[m.group(1).upper(),int(m.group(2))]})
        if OUTPUT_notebook():
            print r'\fn{%s}: sfo1:%0.5f aq:%0.3f swh:%0.3f ns:%d ds: %d rg:%0.1f d1:%0.1f p1:%0.2f pl1:%0.1f'%(filename,v['SFO1'],v['TD']/v['SW_h']/2.0,v['SW_h'],v['NS'],v['DS'],v['RG'],v['D'][1],v['P'][1],v['PL'][1])
            if len(also) > 0:
                for k,val in also.iteritems():
                    if type(val) is list or type(val) is tuple:
                        try:
                            print k,':',lsafe(v[val[0]][val[1]])
                        except:
                            print "(Can't find",k,val,map(type,val),"!)"
                    else:
                        print k,':',lsafe(v[val])
            data = fornotebook.save_data()
            pptvalue = v['SFO1']/data['current_frequency']
            if abs(pptvalue-data['current_ppt'])>1e-4:
                print ('WARNING! ppt value is off!! (',pptvalue,' ppt)')
            else:
                print '%0.4f $ppt$'%pptvalue
            print '\n\n'
#}}}

#{{{regularization
def regularize1d(b,t,tau,alpha):
    """
    What does this get used for?
    Thikinov regularization....
    """
    # for Lx=b
        if size(b) != size(t):
            print "ERROR, size of b doesn't match size of t"
        tau = tau.reshape(1,-1)
        t = t.reshape(-1,1)
        L = exp(-t/tau)
        U,s,V = svd(L,full_matrices=0)
        rms = zeros(len(alpha),dtype='double')
        coeff = zeros((size(tau),size(alpha)),dtype='double')
        fit = zeros((size(t),size(alpha)),dtype='double')
        for j in range(0,len(alpha)):
            S = diag(s / (s**2 + alpha[j]**2))
            x = dot(
                    dot(
                        conj(transpose(V)),
                        dot(S,conj(transpose(U))))
                    ,b)# was b
            fit[:,j] = dot(L,x)
            try:
                coeff[:,j] = x.flatten()
            except:
                print 'shape(coeff)',shape(coeff),'shape(x)',shape(x)
                print 'first',shape(coeff[:,j]),'second',shape(x.reshape(-1,1))
                raise
            rms[j] = linalg.norm(fit[:,j]-b)
        return (coeff,fit,rms)
#}}}


#}}}
