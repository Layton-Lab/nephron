from defs import *
from values import *
import numpy as np
import math


def namgx(cell, memb_id, act, area):

    Km_Nao = 62
    Km_Mgo = 10

    Km_Nai = 40
    Km_Mgi = 0.7

    n_Nao = 2 #1.72
    n_Mgi = 2 #2.4

    m = memb_id[1]

    flux_Nao = cell.conc[0, m] ** n_Nao / (cell.conc[0, m] ** n_Nao + Km_Nao ** n_Nao)

    flux_Nai = Km_Nai / (cell.conc[0, 1] + Km_Nai) #cell.conc[0, 1] / (cell.conc[0, 1] + Km_Nai)

    flux_Mgo = Km_Mgo / (cell.conc[16, m] + Km_Mgo) #cell.conc[16, m] / (cell.conc[16, m] + Km_Mgo)

    flux_Mgi = cell.conc[16, 1] ** n_Mgi / (cell.conc[16, 1] ** n_Mgi + Km_Mgi ** n_Mgi)

    JNaMgx_Mg = area[memb_id[0]][memb_id[1]] * act * flux_Nao * flux_Nai * flux_Mgo * flux_Mgi
    JNaMgx_Na = -2 * JNaMgx_Mg # stoichiometry is Na+:Mg2+ = 2:1

    # print("namgx", JNaMgx_Mg, act, flux_Nao, flux_Nai, flux_Mgo, flux_Mgi)

    return [0, 16], [JNaMgx_Na, JNaMgx_Mg]