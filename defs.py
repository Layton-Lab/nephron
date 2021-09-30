import numpy as np
import collections

solute_id = {
    'Na'    : 0,
    'K'     : 1,
    'Cl'    : 2,
    'HCO3'  : 3,
    'H2CO3' : 4,
    'CO2'   : 5,
    'HPO4'  : 6,
    'H2PO4' : 7,
    'urea'  : 8,
    'NH3'   : 9,
    'NH4'   : 10,
    'H'     : 11,
    'HCO2'  : 12,
    'H2CO2' : 13,
    'glucose' : 14
    }

# valence
solute_val = collections.OrderedDict()
solute_val['Na'] = 1
solute_val['K'] = 1
solute_val['Cl'] = -1
solute_val['HCO3'] = -1
solute_val['H2CO3'] = 0
solute_val['CO2'] = 0
solute_val['HPO4'] = -2
solute_val['H2PO4'] = -1
solute_val['urea'] = 0
solute_val['NH3'] = 0
solute_val['NH4'] = 1
solute_val['H'] = 1
solute_val['HCO2'] = -1
solute_val['H2CO2'] = 0
solute_val['glucose'] = 0

# extract valence values only as an array for convenience
zval = np.asarray(list(solute_val.values()))

compart_id = {
    'Lumen' : 0,
    'Cell'  : 1,
    'ICA'   : 2,
    'ICB'   : 3,
    'LIS'   : 4,
    'Bath'  : 5
    }

# number of compartments
NC = len(compart_id.keys())
# number of solutes
NS = len(solute_id.keys())

class coupled_transport:
    def __init__(self):
        # permeability
        self.perm = 0.0
        # solute id
        self.solute_id = []
        # membrane interface
        self.membrane_id = []
        # coefficient
        self.coef = []

class transporter:
    def __init__(self):
        # activity
        self.act = 0.0
        # which transporter
        self.type = []
        # membrane interface
        self.membrane_id = []

class membrane:
    
    def __init__(self):
        # number of cells
        self.total = 200
        # segment of cell
        self.segment = 'PT'
        # diabete status
        self.diabete = 'Non'
        # unx status
        self.unx = 'N'
        # pregnancy status
        self.preg = 'non'
        # tubule length
        self.len = 0.0
        # luminal diameter
        self.diam = 0
        # surface area
        self.area = np.zeros(NC*NC).reshape((NC,NC))
        # initial surface area
        self.area_init = np.zeros(NC*NC).reshape((NC,NC))
        # reference volume
        self.volref = np.zeros(NC)
        # actual volume
        self.vol = np.zeros(NC)
        # initial volume
        self.vol_init = np.zeros(NC)
        # initial pH
        self.pH = np.zeros(NC)

        # water permeability
        self.dLPV = np.zeros(NC*NC).reshape((NC,NC))
        # solute reflection coefficient
        self.sig = np.ones(NC*NC*NS).reshape((NS,NC,NC))
        # solute permeability
        self.h = np.zeros(NS*NC*NC).reshape((NS,NC,NC))
        # coupled transporters
        self.dLA = []
        # specific transporters
        self.trans = []

        # solute concentrations
        self.conc = np.zeros(NC*NS).reshape((NS,NC))
        # membrane potential
        self.ep = np.zeros(NC)  #[0 for i in range(NC)]

        # reference impermeate concentrations in each compartment
        # denote oncotic pressure in lumen and bath
        self.cimpref = np.zeros(NC)
        # Impermeant Properties
        self.zimp = np.zeros(NC)
        # reference buffer concentrations in each compartment: (concentrations in cell, ICA and ICB are only needed)
        self.cbuftot = np.zeros(NC)


        # pressure
        self.pres = np.zeros(NC)

        # rates used in HCO3/H2CO3 reaction
        self.dkd = np.zeros(NC)
        self.dkh = np.zeros(NC)

        # default sex is male
        self.sex = 'male'

        # default type of nephron is superficial (sup). Another type is juxtamedullary (jux1,jux2,jux3,jux4,jux5)
        self.type = 'sup'

        # inhibition
        self.inhib = ''

        # default is human
        self.species = 'hum'
        self.flag=0
        
        # interstitial concentration parameters of solutes
        # CM: cortico-medullary boundary; OI: outer-inner stripe boundary; Pap: papillary tip
        self.cm = []
        self.oi = []
        self.pap = []
        
        
