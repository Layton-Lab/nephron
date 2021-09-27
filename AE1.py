import numpy as np
from values import *
from defs import *

def AE1(cell,memb_id,xAE1,area):

	bpp = cell.conc[3,2]
	cpp = cell.conc[2,2]
	betapp = bpp/dKbpp
	gampp = cpp/dKcpp

	bp = cell.conc[3,memb_id[1]]
	cp = cell.conc[2,memb_id[1]]
	betap = bp/dKbp
	gamp = cp/dKcp
	xT = xAE1/(1+bpp/172.0)
	summ = (1+betap+gamp)*(Pbpp*betapp+Pcpp*gampp)+(1+betapp+gampp)*(Pbp*betap+Pcp*gamp)
	befflux = cell.area[memb_id[0],memb_id[1]]*xT/summ*(Pbpp*betapp*Pcp*gamp-Pbp*betap*Pcpp*gampp)

	return [3,2],[befflux,-befflux]