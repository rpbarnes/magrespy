
#{{{ k_sigma s(p)
class ksp (fitdata):
    """
    Fit is so slow. do not use this function.
    """
    def __init__(self,*args,**kwargs):
        fitdata.__init__(self,*args,**kwargs)
        #self.symbol_list = [r'ks_{max}',r'p_{half}']
        self.symbol_list = [r'ksmax',r'phalf']
        self.starting_guesses = map(double,[r_[60,0.5],r_[20,0.05],r_[1.,1.]])
        self.guess_lb = r_[0.1,1e-4]
        self.guess_ub = r_[500.,30.]
        self.zero_guess = r_[True,True]
        self.gen_symbolic(r'k_{\sigma} s(p)')
        return
    def fitfunc_raw(self,p,x):
        return p[0]/(p[1]+x)*x
    def fitfunc_raw_symb(self,p,x):
        try:
            return p[0]/(p[1]+x)*x
        except:
            raise CustomError('types of p',p,map(type,p),'type of x',type(x),'x=',x)
#}}}

class emax(fitdata):
    """
    Fit is so slow. do not use this function.
    """
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        newdata = self.copy()
        newdata.sort(self.fit_axis)
        power = newdata.getaxis(self.fit_axis)
        integral = newdata.data
        largest = len(power)
        initial_slope = (integral[largest/4]-integral[0])/(power[largest/4]-power[0])
        approx_emax = integral[-1]/integral[0]
        return [approx_emax,integral[0],-initial_slope] #guess [r'E_{max}',r'v',r'A']
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        #self.function_string = r'$E(p)=v-Apv(1-E_{max})/(1-E_{max}+Ap)$'
        #self.symbol_list = [r'E_{max}',r'v',r'A'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        return (p[1]-(p[2]*x*p[1]*(1.-p[0])/(1.-p[0]+p[2]*x)))
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of E_max, we want 1/(1-E) = 
        '''
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        #print "linfunc passed x=",x,"and y=",y
        rety = 1./(1.-y)
        if yerr != None:
            reterr = yerr/((1.-y)**2)
        mask = isfinite(rety)
        retx = 1./x # for instance, in emax, this is not just x
        xname = r'1 / '+self.fit_axis # same as the fit axis
        yname = r'$\frac{1}{1-E(p)}$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([xname],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        #self.function_string = r'$E(p)=v-Apv(1-E_{max})/(1-E_{max}+Ap)$'
        self.symbol_list = [r'E_{max}',r'v',r'A'] # note that it must not be possible to find part of one of the later strings by searching for one of the earlier strings
        self.gen_symbolic(r'E(p)')
        return

#{{{ filter fitting 
def exp_fit(x,y):
    '''fit an exponential function, disregarding the first point'''
    testpoint = int(len(x)/10)
    startpoint = y[1:6].mean() #because of the filter, sometimes, the very first point can start at 0
    endpoint = y[-testpoint-6:-testpoint].mean()
    initial_slope = (y[testpoint]-startpoint)/(x[testpoint]-x[3])
    # initial slope = -startpoint/t1
    p_ini = [startpoint,-startpoint/initial_slope,endpoint]
    p_out,success = leastsq(exp_errfunc, p_ini,
            args = (x[1:],y[1:]),maxfev=3000,ftol = 1e-4, xtol = 1e-6)
    if (success<0) or (success>4):
        p_out = r_[1.,len(x)*4.,0]
        clf()
        plot(x,exp_fitfunc(p_ini,x))
        plot(x,exp_fitfunc(p_out,x))
        plot(x,y)
        legend(['initial guess','after fit','data'])
        error_plot(success,'There was an error attempting to fit')
    else:
        return p_out

def exp_fitfunc(p,x):
    return p[0]*exp(-x/p[1])+p[2]

def exp_errfunc(p,x,y):
    fit = exp_fitfunc(p,x)
    return fit-y
#}}}
