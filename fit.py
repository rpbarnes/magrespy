
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
