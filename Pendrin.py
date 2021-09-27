import numpy as np
from values import *
from defs import *

def Pendrin(cell,memb_id,xPendrin,area):

	ph=np.zeros(NC)
	for j in range(NC):
		ph[j] = -np.log(abs(cell.conc[11][j])/1.0e3)/np.log(10.0)

	Pbiiclepd = Pbieclepd*Pcliclepd
	Pohiclepd = Poheclepd*Pcliclepd

	ohe = 1.0e-11/cell.conc[11,0]
	ohi = 1.0e-11/cell.conc[11,3]
	alpe = ohe/dKohpd
	alpi = ohi/dKohpd
	game = cell.conc[2,0]/dKclpd
	gami = cell.conc[2,3]/dKclpd
	bete = cell.conc[3,0]/dKbipd
	beti = cell.conc[3,3]/dKbipd

	dele = 1.0+game+bete+alpe
	deli = 1.0+gami+beti+alpi
	etae = 1.0*game+Pbieclepd*bete+Poheclepd*alpe
	etai = Pcliclepd*gami+Pbiiclepd*beti+Pohiclepd*alpi
	sigma = dele*etai+deli*etae

	if cell.pH[3]<7.9:
		adj = np.exp((1/6.4189)*np.log(abs(7.9-ph[3])/0.0115))
	else:
		adj = 1.0

	factpd = float(cell.area[0,3])*xPendrin*Pclepd*adj

	dJclpd = factpd*(1*game*etai-Pcliclepd*gami*etae)/sigma
	dJbipd = factpd*(Pbieclepd*bete*etai-Pbiiclepd*beti*etae)/sigma

	#print(Pbiiclepd,Pohiclepd,ohe,ohi,alpe,alpi,game,gami,bete,beti,dele,deli,etae,etai,sigma,ph[3])
	return [2,3],[dJclpd,dJbipd]#,factpd,cell.area[0,3],xPendrin,Pclepd,adj