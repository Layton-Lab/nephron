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
# newton solvers
#   separate for rat/human, same method with different
#       damping coefficients
#=====================================================
# rat newton solver
def newton_rat(func,x,k,cell):
    if cell.humOrrat != 'rat' and cell.humOrrat != 'mou':
        raise Exception('newton_rat only for rat or mouse model')
    fun=equations.conservation_eqs
    f=np.matrix(fun(x,k))
    TOLpcn = 1
    i = 1
    iter = 0
    maxiter = 150

    # check
    if np.isnan(np.linalg.norm(f)):
        raise Exception('norm(f) is Nan')

    while(np.linalg.norm(f) > 0.0001) and (iter<maxiter+1): #(iter<300)
        if np.linalg.norm(f)>1e12 or np.isnan(np.linalg.norm(f)):
            raise Exception('Newton solver diverged in '+ cell.segment + ' at cell number: ' + str(k))
        elif iter == maxiter:
            print('Warning!!: Newton solver did not converge in <'+str(iter)+' iterations in ' + cell.segment + ' cell number ' + str(k) + '\n')
            print('error size: '+ str(np.linalg.norm(f)))
            print('\n')
            
        i += 1
        J = np.matrix(Jac(fun,x,k))
        IJ = J.I
        F = np.matrix(fun(x,k))
        # PCT
        if cell.segment == 'PT':
            amp = 1.0
        # S3
        elif cell.segment == 'S3':
            amp = 1.0 
        # SDL
        elif cell.segment == 'SDL':
            amp = 1.0
        # LDL
        elif cell.segment == 'LDL':
            if np.linalg.norm(f)>5000:
                amp = 0.5
            else:
                amp = 1.0
        # LAL
        elif cell.segment == 'LAL':
            if np.linalg.norm(f)>5000:
                amp = 0.5
            else:
                amp = 1.0
        # mTAL
        elif cell.segment == 'mTAL':
            if np.linalg.norm(f)>5000:
                amp = 0.5
            elif iter>100:
                amp = 0.95
            else:
                amp = 1.0
        # cTAL
        elif cell.segment == 'cTAL':
            if np.linalg.norm(f)>5000:
                amp = 0.2
            elif np.linalg.norm(f)>1000:
                amp = 0.5
            elif iter>100:
                if np.linalg.norm(f)>1:
                    amp = 0.9
                else:
                    amp = 0.95
            else:
                amp = 1.0
        # DCT
        elif cell.segment == 'DCT':
            amp = 1.0
        # CNT
        elif cell.segment == 'CNT':
            if np.linalg.norm(f)>5000:
                amp = 0.25
            elif np.linalg.norm(f)>1000:
                amp = 0.5
            elif iter > 50:
                if np.linalg.norm(f)>1:
                    amp = 0.8
                else:
                    amp = 0.95
            else:
                amp = 1.0
        # CCD     
        elif cell.segment == 'CCD':
            if np.linalg.norm(f)>5000:
                if k==0:
                    amp = 0.4
                else:
                    amp = 0.5
            elif np.linalg.norm(f)>1000:
                amp = 0.6
            elif iter>75:
                amp = 0.7
            elif iter>50:
                if np.linalg.norm(f)>1:
                    amp = 0.95
                else:
                    amp = 1.0
            else:
                amp = 1.0
        # OMCD
        elif cell.segment == 'OMCD':
            if np.linalg.norm(f)>5000:
                amp = 0.25
            elif np.linalg.norm(f)>1000:
                amp = 0.5
            elif np.linalg.norm(f)>100:
                amp = 1.0
            else:
                amp = 1.0
        # IMCD
        elif cell.segment == 'IMCD':
            if np.linalg.norm(f)>5000:
                if k==0:
                    amp = 0.25 
                else:
                    amp = 0.5 
            elif np.linalg.norm(f)>1000:
                amp = 1.0 #0.5
            elif iter>75:
                if np.linalg.norm(f)>10:
                    amp = 0.75
                elif np.linalg.norm(f)>1:
                    amp = 0.9
                else:
                    amp = 1.0
            elif iter>50:
                if np.linalg.norm(f)>10:
                    amp = 0.8
                else:
                    amp = 0.9
            else:
                amp = 1.0
        else:
            print('What is this segment?', cell.segment)
            raise Exception('cell.segment is not characterized')
        delta = amp*np.array(F*IJ.T)[0]
        x-= delta
        f = np.matrix(fun(x,k))
        iter+=1
        #print(iter, np.linalg.norm(f))
        TOLpcn = np.max(delta/x)
    return x
    


# human newton solver
def newton_human(func,x,k,cell):
    if cell.humOrrat != 'hum':
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
                                amp = 0.13
                            else:
                                amp = 0.17
                        else:
                            if k==0:
                                amp = 0.8
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
        
