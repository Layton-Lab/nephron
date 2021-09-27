from values import *

# across basolateral
def kcc4(conc,memb_id,xKCC4,area):

	bet2 = conc[1][1]/bkkcc
	bet5 = conc[1][4]/bkkcc
	bet6 = conc[1][5]/bkkcc
	gam2 = conc[2][1]/bckcc
	gam5 = conc[2][4]/bckcc
	gam6 = conc[2][5]/bckcc
	dnu2 = conc[10][1]/bmkcc
	dnu5 = conc[10][4]/bmkcc
	dnu6 = conc[10][5]/bmkcc

	sig2 = 1.0+gam2*(1.0+bet2+dnu2)
	sig5 = 1.0+bet5*(1.0+gam5)+dnu5*(1.0+gam5)
	sig6 = 1.0+bet6*(1.0+gam6)+dnu6*(1.0+gam6)

	rho2 = poppkcc+pkccpp*bet2*gam2+pmccpp*dnu2*gam2
	rho5 = popkcc+pkccp*bet5*gam5+pmccp*dnu5*gam5
	rho6 = popkcc+pkccp*bet6*gam6+pmccp*dnu6*gam6

	bigsum5 = sig5*rho2+sig2*rho5
	bigsum6 = sig6*rho2+sig2*rho6

	t1for5 = poppkcc*pkccp*bet5*gam5-popkcc*pkccpp*bet2*gam2
	t2for5 = pmccpp*dnu2*gam2*pkccp*bet5*gam5-pmccp*dnu5*gam5*pkccpp*bet2*gam2
	t3for5 = poppkcc*pmccp*dnu5*gam5-popkcc*pmccpp*dnu2*gam2
	t4for5 = pmccp*dnu5*gam5*pkccpp*bet2*gam2-pmccpp*dnu2*gam2*pkccp*bet5*gam5

	dJk5 = -xKCC4*area[1][4]*(t1for5+t2for5)/bigsum5
	dJm5 = -xKCC4*area[1][4]*(t3for5+t4for5)/bigsum5
	dJc5 = dJk5+dJm5

	t1for6 = poppkcc*pkccp*bet6*gam6-popkcc*pkccpp*bet2*gam2
	t2for6 = pmccpp*dnu2*gam2*pkccp*bet6*gam6-pmccp*dnu6*gam6*pkccpp*bet2*gam2
	t3for6 = poppkcc*pmccp*dnu6*gam6-popkcc*pmccpp*dnu2*gam2
	t4for6 = pmccp*dnu6*gam6*pkccpp*bet2*gam2-pmccpp*dnu2*gam2*pkccp*bet6*gam6

	dJk6 = -xKCC4*area[1][5]*(t1for6+t2for6)/bigsum6
	dJm6 = -xKCC4*area[1][5]*(t3for6+t4for6)/bigsum6
	dJc6 = dJk6 + dJm6

	if memb_id==[1,4]:
		dJk=dJk5
		dJm=dJm5
		dJc=dJc5
	elif memb_id==[1,5]:
		dJk=dJk6
		dJm=dJm6
		dJc=dJc6

	return [1,2,10],[dJk,dJc,dJm]