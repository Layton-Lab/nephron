import equations
import math
import numpy as np
from Newton import Jac

def newton_preg_rat(func,x,k,cell):
    if cell.species != 'rat':
        print('species:' + cell.species)
        raise Exception('newton_preg_rat only for rat model')
    if cell.sex.lower() != 'female':
        print('sex: ' + cell.sex)
        raise Exception('newton_preg only for pregnant female rat')
    fun=equations.conservation_eqs
    f=np.matrix(fun(x,k))
    TOLpcn = 1
    i = 1
    iter = 0
    if k==0:
        maxiter = 300
    else:
        maxiter = 150

    # check
    if np.isnan(np.linalg.norm(f)):
        raise Exception(cell.segment + ' norm(f) is Nan')

    while(np.linalg.norm(f) > 0.0001) and (iter<maxiter+1): 
        if np.linalg.norm(f)>1e15 or np.isnan(np.linalg.norm(f)):
            if cell.segment == 'CCD' or cell.segment == 'OMCD' or cell.segment == 'IMCD':
                raise Exception('Newton solver diverged in ' + cell.segment + ' at cell number: ' + str(k))
            else:
                raise Exception('Newton solver diverged in '+cell.type + ' ' + cell.segment + ' at cell number: ' + str(k))
        elif iter == maxiter:
            if cell.segment == 'CCD' or cell.segment == 'OMCD' or cell.segment == 'IMCD':
                print('Warning!!: Newton solver did not converge in <'+str(iter)+' iterations in ' + cell.segment + ' cell number ' + str(k) + '\n')
            else:
                print('Warning!!: Newton solver did not converge in <'+str(iter)+' iterations in ' + cell.type + ' ' + cell.segment + ' cell number ' + str(k) + '\n')
            print('error size: '+ str(np.linalg.norm(f)) + '\n')

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
            elif iter>75:
                amp = 0.95
            elif iter>100:
                amp = 0.5
            else:
                amp = 1.0
        # LAL
        elif cell.segment == 'LAL':
            if np.linalg.norm(f)>5000:
                amp = 0.5
            elif iter>75:
                amp = 0.95
            else:
                amp = 1.0
        # mTAL
        elif cell.segment == 'mTAL':
            if np.linalg.norm(f)>5000:
                amp = 0.2
            elif np.linalg.norm(f)>2000:
                amp = 0.3
            elif np.linalg.norm(f)>1000:
                amp = 0.5
            elif np.linalg.norm(f)>100:
                amp = 0.9
            else:
                amp = 1.0
        # cTAL
        elif cell.segment == 'cTAL':
            if np.linalg.norm(f)>5000:
                amp = 0.2
            elif np.linalg.norm(f)>2000:
                amp = 0.3
            elif np.linalg.norm(f)>1000:
                amp = 0.5
            elif np.linalg.norm(f)>100:
                amp = 0.9
            else:
                amp = 1.0 
        # DCT
        elif cell.segment == 'DCT':
            # the options were added for the HTN model
            if cell.HT != 'N':
                if np.linalg.norm(f)>1e5:
                    if cell.type == 'jux1':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.4
                    elif cell.type == 'jux4':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.3
                    elif cell.type == 'jux5':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.5
                    elif cell.type == 'jux2':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.4
                    elif cell.type == 'jux3':
                        if k==0:
                            amp=1.0
                        else:
                            if k<160:
                                amp = 0.5 #1.0 #0.25 #0.45 #0.3
                            else:
                                amp = 0.2
                    else:
                        amp = 1.0
                elif np.linalg.norm(f)>5000:
                    if cell.type == 'jux4':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.4 #0.475 #0.5
                    elif cell.type == 'jux1':
                        if k==0:
                            amp = 1.0
                        elif k>180:
                            amp = 0.3 #1.0
                        else:
                            amp = 0.25 #0.45 #0.3
                    elif cell.type == 'jux2':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.3
                    elif cell.type == 'jux3':
                        if k==0:
                            amp = 1.0
                        else:
                            if k<160:
                                amp = 0.75 #0.5 #1.0 #0.225 #0.3 #0.4
                            else:
                                amp = 0.25
                    else:
                        amp = 1.0
                elif np.linalg.norm(f)>1000:
                    if cell.type == 'jux1':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.4
                    elif cell.type == 'jux2':
                        if k==0:
                            amp = 1.0
                        elif k>160:
                            amp = 0.5
                        else:
                            amp = 1.0
                    elif cell.type == 'jux3':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.75 #0.4
                    elif cell.type == 'jux4':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.5
                    else:
                        amp = 1.0
                elif np.linalg.norm(f)>250:
                    if cell.type == 'jux1':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 1.0
                    elif cell.type == 'jux4':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.7
                    elif cell.type == 'jux3':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 0.7
                    else:
                        amp = 1.0
                elif np.linalg.norm(f)>50:
                    if cell.type == 'jux4':
                        if k==0:
                            amp = 1.0
                        else:
                            amp = 1.0 #0.7 #0.8
                    elif cell.type == 'jux3':
                        if k==0:
                            amp = 1.0
                        elif k<160:
                            amp = 0.7
                        else:
                            amp = 0.25 #0.5 #1.0
                    else:
                        amp = 1.0
                elif iter>20:
                    if cell.type == 'jux4':
                        if np.linalg.norm(f)>10:
                            amp = 0.75
                        else:
                            amp = 1.0
                    elif cell.type == 'jux1':
                        if np.linalg.norm(f)>10:
                            amp = 0.75
                        else:
                            amp = 1.0
                    else:
                        amp = 1.0
                else:
                    amp = 1.0
            else:
                amp = 1.0
        # CNT
        elif cell.segment == 'CNT':
            if np.linalg.norm(f)>1e6:
                if k==0:
                    amp = 0.4
                else:
                    amp = 0.35 #0.55
            elif np.linalg.norm(f)>5000:
                if cell.preg == 'mid':
                    if cell.type == 'jux1':
                        if cell.HT != 'N':
                            if k==0:
                                amp = 1.0
                            else:
                                amp = 0.3
                        else:
                            amp = 0.275
                    elif cell.type == 'jux2':
                        if cell.HT != 'N':
                            amp = 0.25
                        else:
                            amp = 0.3
                    elif cell.type == 'jux3':
                        if k==0:
                            amp = 0.5
                        else:
                            amp = 0.25
                    elif cell.type == 'jux4':
                        amp = 0.25 #0.277
                    elif cell.type == 'jux5':
                        if k==0:
                            amp = 1.0 #0.4
                        else:
                            amp = 0.25 #0.5
                    else:
                        amp = 0.5
                elif cell.preg == 'late':
                    if cell.type == 'jux1':
                        amp = 0.45
                    elif cell.type == 'jux2':
                        amp = 0.4
                    elif cell.type == 'jux3':
                        amp = 0.4
                    elif cell.type == 'jux4':
                        amp = 0.5
                    else:
                        amp = 0.65
            elif np.linalg.norm(f)>2000:
                if cell.preg == 'mid':
                    amp = 0.85
                elif cell.preg == 'late':
                    amp = 0.85
            elif np.linalg.norm(f)>1000:
                if cell.preg == 'mid':
                    amp = 0.9 #1.0 #0.8
                elif cell.preg == 'late':
                    amp = 0.9
            elif np.linalg.norm(f)>150:
                if cell.preg == 'mid':
                    amp = 1.0
                elif cell.preg == 'late':
                    amp = 0.95
            elif iter>100:
                amp = 0.55
            elif iter > 30:
                if np.linalg.norm(f)>15:
                    if cell.type == 'jux4':
                        amp = 1.0
                    else:
                        amp = 0.85
                elif np.linalg.norm(f)>5:
                    if cell.type == 'jux4':
                        amp = 0.75
                    else:
                        amp = 0.5
                else:
                    if cell.type == 'jux4':
                        amp = 0.5
                    else:
                        amp = 0.76
            elif iter>10:
                if np.linalg.norm(f)<1:
                    amp = 0.5
                else:
                    amp = 1.0
            else:
                amp = 1.0 
        # CCD     
        elif cell.segment == 'CCD':
            if np.linalg.norm(f)>1e6:
                amp = 0.3
            elif cell.preg == 'mid':
                if np.linalg.norm(f)>5000:
                    if k==0:
                        if cell.HT == 'Y':
                            amp = 1.0
                        else:
                            amp = 0.5
                    else:
                        amp = 0.75
                elif np.linalg.norm(f)>1000:
                    amp = 0.75   
                elif iter > 50:
                    if np.linalg.norm(f)>1:
                        amp = 0.8
                    else:
                        amp = 0.95
                else:
                    amp = 1.0
            elif cell.preg == 'late':
                if np.linalg.norm(f)>5000:
                    if k==0:
                        amp = 0.35
                    else:
                        amp = 0.55
                elif np.linalg.norm(f)>2000:
                    if k==0:
                        amp = 0.5
                    else:
                        amp = 0.75 
                elif np.linalg.norm(f)>1000:
                    if k==0:
                        amp = 0.75
                    else:
                        amp = 0.95
                elif iter > 100:
                    amp = 1.0
                elif iter > 125:
                    if np.linalg.norm(f)>1:
                        amp = 0.8
                    else:
                        amp = 0.75
                else:
                    amp = 1.0
        # OMCD
        elif cell.segment == 'OMCD':
            if np.linalg.norm(f)>1e6:
                amp = 0.4
            elif np.linalg.norm(f)>5000:
                if k==0:
                    amp = 0.75
                else:
                    amp = 1.0
            elif np.linalg.norm(f)>1000:
                if k==0:
                    amp = 0.75
                else:
                    amp = 1.0
            elif np.linalg.norm(f)>100:
                amp = 1.0
            elif iter>50:
                amp = 0.9
            elif iter>100:
                amp = 0.8
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
            elif iter>100:
                if np.linalg.norm(f)>10:
                    amp = 1.0
                elif np.linalg.norm(f)>1:
                    amp = 0.8
                else:
                    amp = 0.5
            elif iter>75:
                if np.linalg.norm(f)>10:
                    amp = 1.0
                elif np.linalg.norm(f)>5:
                    amp = 0.85
                elif np.linalg.norm(f)>1:
                    amp = 0.62
                else:
                    amp = 0.5
            elif iter>50:
                if np.linalg.norm(f)>10:
                    amp = 0.8
                elif np.linalg.norm(f)>1:
                    amp = 0.75
                else:
                    amp = 1.0
            else:
                amp = 1.0
        else:
            print('What is this segment?', cell.segment)
            raise Exception('cell.segment is not characterized')
        delta = amp*np.array(F*IJ.T)[0]
        xold = x
        x -= delta
        f = np.matrix(fun(x,k))

        # extra check
        if np.isnan(np.linalg.norm(f)) or np.linalg.norm(f)>1e14:
            amp = amp/2 
            delta = amp*np.array(F*IJ.T)[0]
            x = xold - delta
            f = np.matrix(fun(x,k))

        iter+=1
        #print(iter, np.linalg.norm(f))
        TOLpcn = np.max(delta/x)
        #print(f)
        #input('pausing...')
    return x
        
