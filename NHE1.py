import numpy as np
from values import *
from defs import *

def NHE1(cell,memb_id,NHE1,area):

	affnao = 34.0
	affnai = 102.0
	affho = 0.0183e-3
	affhi = 0.054e-3

	nai = cell.conc[0,memb_id[0]]
	hi = cell.conc[11,memb_id[0]]
	nao = cell.conc[0,memb_id[1]]
	ho = cell.conc[11,memb_id[1]]

	Fno = (nao/affnao)/(1+nao/affnao+ho/affho)
	Fni = (nai/affnai)/(1+nai/affnai+hi/affhi)
	Fho = (ho/affho)/(1+nao/affnao+ho/affho)
	Fhi = (hi/affhi)/(1+nai/affnai+hi/affhi)

	E2mod1 = (Fni+Fhi)/(Fni+Fhi+Fno+Fho)
	E1mod1 = 1-E2mod1
	E2mod2 = (Fni**2+Fhi**2)/(Fni**2+Fhi**2+Fno**2+Fho**2)
	E1mod2 = 1-E2mod2
	Fmod1 = (hi**2)/(hi**2+(0.3e-3)**2)
	Rnhe = 1.0e3*(1*(1-Fmod1)*(E2mod2*(Fno**2)-E1mod2*(Fni**2))+Fmod1*(E2mod1*Fno-E1mod1*Fni))

	dJNaH = -area[memb_id[0],memb_id[1]]*NHE1*Rnhe
	#print(Fmod1,E2mod2,Fno,E1mod2,Fni,E2mod1,E1mod1)
	#print(nai,hi,nao,ho)
	return [0,11],[dJNaH,-dJNaH]