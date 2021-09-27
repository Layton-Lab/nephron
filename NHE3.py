# NHE3 exchanger modules
#---------------------------------------------------------------------
# NHE3 EXCHANGER AT LUMINAL MEMBRANE OF CELL
# n is for Na, c is for Cl
# p (or prime) is for the luminal compartment (M)
# pp (or double prime) is for the cytosolic compartment (I)
#---------------------------------------------------------------------
#   conc[NS][NC]: all concentrations in all compartments
#   ep[NC]:       all potential in all compartments
#   memb_id:      specifies interface, (lumen,cell), (cell,bath), or (cell,LIS)
#   CT:           activity level
#   area:         surface area of all membranes
#
#   returns:      solute IDs and corresponding fluxes
#

# parameters
PaNH = 8000.0
PbNH = 8000.0*0.48/1.6
PcNH = 8000.0
dKaNH = 30.0
dKbNH = 72.0e-6
dKcNH = 27.0
fMNH = 2.0
dKINH = 1.0e-3

def nhe3(cell,ep,memb_id,xNHE3,area):
	
	ap = cell.conc[0][memb_id[0]]  # Na+
	bp = cell.conc[11][memb_id[0]]  # H+
	cp = cell.conc[10][memb_id[0]]  # NH4+
	app = cell.conc[0][memb_id[1]]
	bpp = cell.conc[11][memb_id[1]]
	cpp = cell.conc[10][memb_id[1]]

	alp = ap / dKaNH
	alpp = app / dKaNH
	betap = bp / dKbNH
	betapp = bpp / dKbNH
	gamp = cp / dKcNH
	gampp = cpp / dKcNH
	fmod = fMNH*cell.conc[11][memb_id[1]]/(cell.conc[11][memb_id[1]]+dKINH)
	sum1 = (1+alp+betap+gamp)*(PaNH*alpp+PbNH*betapp+PcNH*gampp)
	sum2 = (1+alpp+betapp+gampp)*(PaNH*alp+PbNH*betap+PcNH*gamp)
	sum = sum1 + sum2
	termNaH = fmod*PaNH*PbNH*(alp*betapp-alpp*betap)
	termNaNH4 = fmod*PaNH*PcNH*(alp*gampp-alpp*gamp)
	termHNH4 = fmod*PbNH*PcNH*(betap*gampp-betapp*gamp)
	dJNHEsod = area[memb_id[0]][memb_id[1]]*xNHE3*(termNaH+termNaNH4)/sum
	dJNHEprot = area[memb_id[0]][memb_id[1]]*xNHE3*(-termNaH+termHNH4)/sum
	dJNHEamm = area[memb_id[0]][memb_id[1]]*xNHE3*(-termNaNH4-termHNH4)/sum

	if cell.segment == 'PT' or cell.segment == 'S3':
	# The rate constants are different for the PT. They need to be divided by 8000./792.
		dJNHEsod = dJNHEsod*792.0/8000.0
		dJNHEprot = dJNHEprot*792.0/8000.0
		dJNHEamm = dJNHEamm*792.0/8000.0
	   
	return [0,11,10],[dJNHEsod,dJNHEprot,dJNHEamm]
