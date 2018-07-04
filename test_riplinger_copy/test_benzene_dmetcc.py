import scipy
import numpy,os
from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb

from pyscf import gto,scf,cc,mp, symm, ao2mo, mp

atoms=[
['C',(0.000, 1.396, 0.000)],
['C',(1.209, 0.698, 0.000)],
['C',(1.209, -0.698,0.000)],
['C',(0.000, -1.396,0.000)],
['C',(-1.209, -0.698,0.000)],
['C',(-1.209, 0.698,0.000)],
['H',(0.000, 2.479,0.000)],
['H',(2.147, 1.240,0.000)],
['H',(2.147, -1.240,0.000)],
['H',(0.000, -2.479,0.000)],
['H',(-2.147, -1.240,0.000)],
['H',(-2.147, 1.240,0.000)]
]


mol = gto.M(atom = atoms, basis='cc-pvdz') # read the geometry, in A
#m = scf.RHF(mol)
#m.kernel()

#mycc = cc.CCSD(m)
#mycc.kernel()

#mymp2 = mp.MP2(m)
#mymp2.kernel()
del mol

bs     = 'dz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments = [[0],[6],[1],[7],[2],[8],[3],[9],[4],[5],[11],[10]]
fragment_spins = [0,0,0,0,0,0,0,0,0,0,0,0]

thresh   = 1.0e-6
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'cc'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)

