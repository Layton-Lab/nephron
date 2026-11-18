#from numba import jit
import equations
import math
import numpy as np

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def Jac(func,x,k):

    epsfcn = 1.0e-3
    epsmch = 1.0e-3
    eps = math.sqrt(max(epsfcn,epsmch))
    
    Jfun=[[0 for i in range(len(x))] for i in range (len(x))]
    
    
    wa1=func(x,k)
    for i in range(len(x)):
        temp=x[i]
        h=eps*abs(temp)
        if (h==0):
            h=eps
        x[i]=temp+h
        fvec=func(x,k)
        x[i]=temp
        for j in range(len(x)):
            Jfun[j][i]=(-wa1[j]+fvec[j])/h
    
    return Jfun
    
#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit       
def broyden(func,x,k,type):
    fun=equations.conservation_eqs
    f=np.matrix(fun(x,k))
    J=np.matrix(Jac(fun,x,k))
#    IJ=np.linalg.inv(J)
    IJ=J.I
    dx=np.ones(x.shape)
    i=0
    iter=0
    #while(np.max(dx)>0.0001):
    while(np.linalg.norm(f)>0.0001) and (iter<500):
        
        f_previous=f
        x_previous=x
        
        x=x-np.array(f*IJ.T)[0]
        
        f=np.matrix(fun(x,k))
#                print('x',x)
#                print('f',f)

        df=f-f_previous
        dx=x-x_previous
        # #
        # #-------------------------------------------------------
        # #using good broyden
        # dx=np.array([dx])
        # df=np.array([df])
        # dx=dx.T
        # df=df.T

        # IJ = IJ+(dx-IJ*df)*(dx.T*IJ)/np.inner(dx.T*IJ,df)
        # #-------------------------------------------------------

        IJ=IJ-np.outer(IJ*f.T,dx)*IJ/np.inner(dx,dx+(IJ*f.T).T)
        #print(i)
        iter+=1
        #print(iter,np.linalg.norm(f))
        
#        input("Pausing! Press Enter to continue...")
        #J=J+np.outer((df-dx*J.T),dx)/np.linalg.norm(x)**2
        #IJ=J.I

    return x

#=====================================================
# Newton solver for the rat and mouse models (the human solver below is separate)
#=====================================================
# One Newton step is  x -= amp * J^-1 F  with a numerical Jacobian.  The damping factor amp is a hand-tuned function of the segment, the
# cell index k, the iteration count and the residual norm (see _damping).  On top of that a positivity safeguard (fraction-to-boundary
# rule) limits a step that would drive a concentration below POSFRAC of its current value; it is switched on per segment:
#   POS_FIRST_CELL   first cell (k == 0) only.  S3: needed for a halved PT Na-K-ATPase; CCD: needed for the first cell.
#   POS_EVERY_CELL   every cell.  CCD: needed for the high-flow diabetic mouse cases (SGLT2 inhibition + NHE3 down), where a full step
#                    drove a ~3e-5 mM concentration negative at cell 1.  NOT used for the OMCD (its valid solution has H2CO2 ~ -1e-4).
#   POS_IF_NEG_FORMATE  first cell is solved plainly and re-solved WITH the safeguard only if the solution has a negative formate
#                    (some female-rat OMCD cases otherwise stay on the negative-formate branch).
# Environment switches:  NEPH_POSSAFE=0 turns every safeguard off (original behaviour);  NEPH_STRICT=1 raises instead of continuing
# when the iteration cap is reached.  Results of every segment/species without a safeguard are identical to the pre-safeguard code.
import os, time as _time
STATS = []        # one record per solve: dict(seg, k, iters, resid, sec); read by the run scripts to count non-converged cells
STRICT = os.environ.get('NEPH_STRICT','0') == '1'
POSSAFE = os.environ.get('NEPH_POSSAFE','1') != '0'
POSFRAC = 0.1
POS_FIRST_CELL = {'CCD', 'S3'}
POS_EVERY_CELL = {'CCD'}
POS_IF_NEG_FORMATE = {'OMCD'}
CD_SEGMENTS = ('CCD', 'OMCD', 'IMCD')   # segments whose error messages do not name the nephron type

class NewtonNotConverged(Exception):
    pass

def _where(cell):
    return cell.segment if cell.segment in CD_SEGMENTS else cell.type + ' ' + cell.segment

def _damping(seg, k, it, nf):
    """Damping factor of the Newton step: segment, cell index k, iteration count it, residual norm nf.
    Segments not listed (PT, S3, SDL, LAL, DCT) take the full step."""
    if seg in ('PT', 'S3', 'SDL', 'LAL', 'DCT'):
        return 1.0
    if seg == 'LDL':
        return 0.9 if it > 75 else 1.0
    if seg == 'mTAL':
        if nf > 1500:
            return 0.2 if k == 0 else 0.5
        if nf > 100:
            return 0.2 if k == 0 else 0.95
        return 0.95 if it > 100 else 1.0
    if seg == 'cTAL':
        if nf > 5000: return 0.2
        if nf > 2000: return 0.3
        if nf > 1000: return 0.5
        if nf > 100:  return 0.9
        if it > 100:  return 0.9 if nf > 1 else 0.95
        return 1.0
    if seg == 'CNT':
        return 0.3
    if seg == 'CCD':
        if nf > 1e8:  return 0.25
        if nf > 5000: return 0.5 if k == 0 else 0.8
        if nf > 1000: return 0.7 if k == 0 else 0.8
        if it > 75:   return 0.75
        if it > 50:   return 0.8 if nf > 1 else 0.95
        if it > 30:   return 0.8 if nf > 20 else (0.9 if nf > 5 else 1.0)
        return 1.0
    if seg == 'OMCD':
        if k == 0:    return 0.4
        if nf > 1e6:  return 0.45
        if nf > 1e4:  return 0.65
        if nf > 100:  return 1.0
        if it > 115:  return 0.5
        if it > 75:   return 0.65
        return 0.7
    if seg == 'IMCD':
        if k == 0:
            return 0.12
        return 1.0 if nf < 100 else 0.5   # undamped once close to the solution
    raise Exception('cell.segment is not characterized: ' + seg)

def _limit_step(x, delta):
    """Fraction-to-boundary: scale the step so that no concentration (the leading entries of x: 5*NS for CNT/CCD/OMCD, 3*NS for the
    three-compartment segments) falls below POSFRAC of its current value."""
    nconc = 75 if len(x) == 86 else 45
    xc = np.asarray(x, dtype=float)[:nconc]; dc = delta[:nconc]
    mask = (dc > 0) & (xc > 0)
    if mask.any():
        smax = np.min((1.0-POSFRAC)*xc[mask]/dc[mask])
        if smax < 1.0:
            delta = delta*max(smax, 1e-3)
    return delta

def _newton_rat(func, x, k, cell, force_safeguard=False):
    if cell.species != 'rat' and cell.species != 'mou':
        raise Exception('newton_rat only for rat or mouse model')
    _t0 = _time.time()
    fun = equations.conservation_eqs
    f = np.matrix(fun(x,k))
    iter = 0
    maxiter = 300 if k == 0 else 150
    safeguard = POSSAFE and ((k == 0 and (cell.segment in POS_FIRST_CELL or force_safeguard)) or cell.segment in POS_EVERY_CELL)

    if np.isnan(np.linalg.norm(f)):
        raise Exception(cell.segment + ' norm(f) is Nan' if cell.segment in CD_SEGMENTS else cell.segment + ' ' + cell.type + ' norm(f) is Nan')

    while (np.linalg.norm(f) > 0.0001) and (iter < maxiter+1):
        nf = np.linalg.norm(f)
        if nf > 1e14 or np.isnan(nf):
            raise Exception('Newton solver diverged in ' + _where(cell) + ' at cell number: ' + str(k))
        elif iter == maxiter:
            print('Warning!!: Newton solver did not converge in <'+str(iter)+' iterations in ' + _where(cell) + ' cell number ' + str(k) + '\n')
            print('error size: '+ str(nf))
            print('\n')
            if STRICT:
                raise NewtonNotConverged(cell.segment + ' cell ' + str(k) + ' |f|=' + str(nf))

        J = np.matrix(Jac(fun,x,k))
        IJ = J.I
        F = np.matrix(fun(x,k))     # re-evaluated on purpose: conservation_eqs keeps state in the cell objects, and Jac leaves it at the last perturbation
        amp = _damping(cell.segment, k, iter, nf)
        delta = amp*np.array(F*IJ.T)[0]
        if safeguard:
            delta = _limit_step(x, delta)
        x -= delta
        f = np.matrix(fun(x,k))
        iter += 1
    STATS.append(dict(seg=cell.segment, k=int(k), iters=int(iter), resid=float(np.linalg.norm(f)), sec=_time.time()-_t0))
    return x

def newton_rat(func,x,k,cell):
    """First-cell wrapper for the segments of POS_IF_NEG_FORMATE: solve cell 0 plainly and, only if any formate concentration
    (solute 12) is negative in the solution, re-solve it from the same start with the positivity safeguard."""
    if POSSAFE and k == 0 and cell.segment in POS_IF_NEG_FORMATE and cell.segment not in POS_FIRST_CELL:
        x0 = np.array(x, copy=True)
        sol = _newton_rat(func,x,k,cell)
        # formate is solute 12; 5-compartment layout (CNT/CCD/OMCD) is solute-major x[5*j+c], 3-compartment is compartment-major x[15*c+j]
        sol_a = np.array(sol, dtype=float)
        form = sol_a[60:65] if len(sol_a) > 60 else sol_a[[12,27,42]]
        if form.min() < 0:
            print('Note: negative formate in the first '+cell.segment+' cell (min %.3g); re-solving it with the positivity safeguard\n' % form.min())
            sol = _newton_rat(func,x0,k,cell,force_safeguard=True)
        return sol
    return _newton_rat(func,x,k,cell)

# human newton solver
def newton_human(func,x,k,cell):
    if cell.species != 'hum':
        raise Exception('newton_human only for human model')
    fun=equations.conservation_eqs
    f=np.matrix(fun(x,k))
    TOLpcn = 1
    i = 1
    iter = 0
    while(np.linalg.norm(f) > 0.0001) and (iter<150): #(iter<300)
        i += 1
        J = np.matrix(Jac(fun,x,k))
        IJ = J.I
        F = np.matrix(fun(x,k))
        # PCT
        if cell.segment == 'PT':
            amp = 1.0
        # S3
        elif cell.segment == 'S3':
            if cell.diabete == 'Non' and cell.unx == 'Y':
                if cell.sex == 'male':
                    amp = 1.0
                if cell.sex == 'female':
                    amp = 0.6
            else:
                amp = 1.0
        # SDL
        elif cell.segment == 'SDL':
            amp = 1.0
        # LDL
        elif cell.segment == 'LDL':
            amp = 1.0
        # LAL
        elif cell.segment == 'LAL':
            amp = 1.0
        # mTAL
        elif cell.segment == 'mTAL':
            if cell.inhib == 'ACE':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 0.8
            elif cell.diabete == 'Non':
                if np.linalg.norm(f)>100:
                    amp=0.2
                else:
                    amp=1.0
            elif cell.diabete == 'Moderate':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            elif cell.diabete == 'Severe':
                if np.linalg.norm(f)>100:
                    amp = 0.2
                else:
                    amp = 1.0
            else:
                amp = 1.0
        elif cell.segment == 'cTAL':
            if np.linalg.norm(f)>100:
                amp = 0.2
            else:
                amp = 1.0
        # MD
        elif cell.segment == 'MD':
            if np.linalg.norm(f)>100:
                amp = 0.2
            else:
                amp = 0.8
        # DCT
        elif cell.segment == 'DCT':
            if cell.inhib == 'ACE':
                if cell.sex == 'male':
                    amp = 0.5
                elif cell.sex == 'female':
                    amp = 0.8
            else:
                if cell.sex == 'female':
                    if cell.type == 'jux1':
                        amp = 0.9
                    elif cell.type == 'jux2':
                        amp = 0.7
                    elif cell.type == 'jux3':
                        amp = 0.7
                    elif cell.type == 'sup':
                        if np.linalg.norm(f)>2000:
                            amp = 0.5
                        else:
                            amp = 0.5
                    else:
                        amp = 1.0
                else:
                    amp = 0.5      
        # CNT
        elif cell.segment == 'CNT':
            if np.linalg.norm(f)>5000:
                if k==0:
                    amp = 0.2
                else:
                    amp = 0.3
            elif np.linalg.norm(f)>1000:
                if k==0:
                    amp = 0.5
                else:
                    amp = 0.7
            else:
                amp = 1.0
        # CCD
        elif cell.segment == 'CCD':
            if np.linalg.norm(f)>5000:
                if k==0:
                    amp = 0.25
                else:
                    amp = 0.5
            elif np.linalg.norm(f)>1000:
                amp = 0.75
            elif iter>50:
                if np.linalg.norm(f)>1:
                    amp = 0.8
                else: 
                    amp = 0.95
            else:
                amp = 1.0

        # OMCD
        elif cell.segment == 'OMCD':
            if np.linalg.norm(f)>100:
                amp = 0.5
            else:
                amp = 0.8
        # IMCD
        elif cell.segment == 'IMCD':
            if cell.inhib == 'ACE':
                if np.linalg.norm(f)>100:
                    if k==0:
                        amp=0.2
                    else:
                        amp = 0.1
                else:
                    if k==0:
                        amp = 1.5
                    else:
                        amp = 0.9
            elif cell.diabete == 'Non' and cell.unx == 'N':
                if cell.sex == 'female':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.2
                        else:
                            amp = 0.1
                    else:
                        if k==0:
                            amp = 0.5
                        else:
                            amp = 0.5
                elif cell.sex == 'male':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.27
                        else:
                            amp = 0.2
                    else:
                        if k==0:
                            amp = 0.8
                        else:
                            amp = 0.8
            elif cell.diabete == 'Non' and cell.inhib == 'SGLT2' and cell.unx == 'Y':
                if cell.sex == 'female':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.2
                        else:
                            amp = 0.1
                    else:
                        if k==0:
                            amp = 0.5
                        else:
                            amp = 0.5
                elif cell.sex == 'male':
                    if np.linalg.norm(f)>100:
                        if k==0:
                            amp = 0.23
                        else:
                            amp = 0.8
                    else:
                        amp = 0.8
            elif cell.diabete == 'Moderate':
                if cell.inhib != 'SGLT2':
                    if cell.sex == 'female':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.2
                            else:
                                amp = 0.1
                        else:
                            amp = 0.5
                    elif cell.sex == 'male':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.19
                            else:
                                amp = 0.17
                        else:
                            if k==0:
                                amp = 0.8
                            else:
                                amp = 0.8
                elif cell.inhib == 'SGLT2':
                    if cell.sex == 'male':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.2
                        else:
                            amp = 0.8
                    elif cell.sex == 'female':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.7
                            else:
                                amp = 0.2
                        else:
                            amp = 0.8
            elif cell.diabete == 'Severe':
                if cell.inhib == 'SGLT2':
                    if cell.sex == 'female':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.2
                            else:
                                amp = 0.1
                        else:
                            if k==0:
                                amp = 0.5
                            else:
                                amp = 0.5
                    elif cell.sex == 'male':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.23
                            else:
                                amp = 0.17
                        else:
                            if k==0:
                                amp = 0.8
                            else:
                                amp = 0.8
                else:
                    if cell.sex == 'female':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.2
                            else:
                                amp = 0.1
                        else:
                            if k == 0:
                                amp = 0.5
                            else:
                                amp = 0.5
                    elif cell.sex == 'male':
                        if np.linalg.norm(f)>100:
                            if k==0:
                                amp = 0.1
                            else:
                                amp = 0.17
                        else:
                            if k==0:
                                amp = 0.2
                            else:
                                amp = 0.8
            else:
                amp = 1.0
        else:
            print('What is this segment?', cell.segment)
            raise Exception('cell.segment is not characterized')
        delta = amp*np.array(F*IJ.T)[0]
        x -= delta
        f = np.matrix(fun(x,k))
        iter += 1
        #print(iter, np.linalg.norm(f))
        TOLpcn = np.max(delta / x)
    return x
        
