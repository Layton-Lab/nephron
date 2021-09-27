from defs import *
from values import *
import numpy as np

def compute_cotransport (cell,delmu,jsol):

    numLA = len(cell.dLA)
    flux = np.zeros(numLA)
    for i in range(numLA):
        
        sid = list(cell.dLA[i].solute_id)
        mid = cell.dLA[i].membrane_id

        this_delmu = list(delmu[i][mid[0]][mid[1]] for i in sid)
        flux = cell.area[mid[0]][mid[1]] * cell.dLA[i].perm * sum(np.array(cell.dLA[i].coef) * this_delmu)

        ind = 0
        
        for k in sid:
            jsol[k][mid[0]][mid[1]] += cell.dLA[i].coef[ind]*flux         
            ind += 1        

    return jsol
