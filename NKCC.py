from values import *
import numpy as np
from defs import *

def nkcc2(cell,memb_id,xNKCC2,area,isoform):
	if cell.diabete != 'Non':
		if cell.segment == 'mTAL':
			xNKCC2 = xNKCC2*1.1

	bn2_dic = {'F':bn2F,'A':bn2A,'B':bn2B}
	bc2_dic = {'F':bc2F,'A':bc2A,'B':bc2B}
	bk2_dic = {'F':bk2F,'A':bk2A,'B':bk2B}
	bm2_dic = {'F':bm2F,'A':bm2A,'B':bm2B}
	poppnkcc_dic = {'F':poppnkccF,'A':poppnkccA,'B':poppnkccB}
	popnkcc_dic = {'F':popnkccF,'A':popnkccA,'B':popnkccB}
	pnkccp_dic = {'F':pnkccpF,'A':pnkccpA,'B':pnkccpB}
	pnkccpp_dic = {'F':pnkccppF,'A':pnkccppA,'B':pnkccppB}
	pnmccp_dic = {'F':pnmccpF,'A':pnmccpA,'B':pnmccpB}
	pnmccpp_dic = {'F':pnmccppF,'A':pnmccppA,'B':pnmccppB}

	bn2 = bn2_dic[isoform]
	bc2 = bc2_dic[isoform]
	bk2 = bk2_dic[isoform]
	bm2 = bm2_dic[isoform]
	poppnkcc = poppnkcc_dic[isoform]
	popnkcc = popnkcc_dic[isoform]
	pnkccp = pnkccp_dic[isoform]
	pnkccpp = pnkccpp_dic[isoform]
	pnmccp = pnmccp_dic[isoform]
	pnmccpp = pnmccpp_dic[isoform]

	alp1 = cell.conc[0][memb_id[0]]/bn2
	alp2 = cell.conc[0][memb_id[1]]/bn2
	bet1 = cell.conc[1][memb_id[0]]/bk2
	bet2 = cell.conc[1][memb_id[1]]/bk2
	gam1 = cell.conc[2][memb_id[0]]/bc2
	gam2 = cell.conc[2][memb_id[1]]/bc2
	dnu1 = cell.conc[10][memb_id[0]]/bm2
	dnu2 = cell.conc[10][memb_id[1]]/bm2

	sig1=1+alp1+alp1*gam1*(1+bet1+bet1*gam1+dnu1+dnu1*gam1)
	sig2=1+gam2*(1+bet2+bet2*gam2+bet2*gam2*alp2+dnu2+dnu2*gam2+dnu2*gam2*alp2)
	rho1=popnkcc+pnkccp*alp1*bet1*(gam1**2)+pnmccp*alp1*dnu1*(gam1**2)
	rho2=poppnkcc+pnkccpp*alp2*bet2*(gam2**2)+pnmccpp*alp2*dnu2*(gam2**2)

	bigsum = sig1*rho2+sig2*rho1

	t1 = poppnkcc*pnkccp*alp1*bet1*(gam1**2)-popnkcc*pnkccpp*alp2*bet2*(gam2**2)
	t2 = poppnkcc*pnmccp*alp1*dnu1*(gam1**2)-popnkcc*pnmccpp*alp2*dnu2*(gam2**2)
	t3 = pnmccpp*alp2*dnu2*(gam2**2)*pnkccp*alp1*bet1*(gam1**2)
	t4 = pnmccp*alp1*dnu1*(gam1**2)*pnkccpp*alp2*bet2*(gam2**2)

	dJnNKCC2 = xNKCC2*area[memb_id[0]][memb_id[1]]*(t1+t2)/bigsum
	dJkNKCC2 = xNKCC2*area[memb_id[0]][memb_id[1]]*(t1+t3-t4)/bigsum
	dJmNKCC2 = xNKCC2*area[memb_id[0]][memb_id[1]]*(t2+t4-t3)/bigsum
	dJcNKCC2 = 2*dJnNKCC2

	return [0,1,2,10],[dJnNKCC2,dJkNKCC2,dJcNKCC2,dJmNKCC2]

def nkcc1(cell,memb_id,xNaKCl2,delmu):

	dJNaKCl2 = cell.area[memb_id[0],memb_id[1]]*xNaKCl2*(delmu[0,memb_id[0],memb_id[1]]+delmu[1,memb_id[0],memb_id[1]]+2*delmu[2,memb_id[0],memb_id[1]])

	return [0,1,2],[dJNaKCl2,dJNaKCl2,2*dJNaKCl2]