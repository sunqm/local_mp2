import scipy
import numpy,os
from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb

from pyscf import gto,scf,cc,mp, symm, ao2mo, mp

atoms=[
['C',(-0.64866, 1.22109, 1.54261)],
['C',(1.41161, 0.05330, 1.39459)],
['C',(-0.74501, -1.17043, 1.54726)],
['C',(0.64866, -1.22109, 1.54261)],
['C',(1.41161, -0.05330, 1.39459)],
['C',(2.79287, -0.11965, 0.79362)],
['C',(-2.79287, 0.11965, 0.79362)],
['C',(2.79287, 0.11965, -0.79362)],
['C',(-2.79287, -0.11965, -0.79362)],
['C',(1.41161, 0.05330, -1.39459)],
['C',(0.74501, -1.17043, -1.54726)],
['C',(-0.64866, -1.22109, -1.54261)],
['C',(-1.41161, -0.05330, -1.39459)],
['C',(-0.74501, 1.17043, -1.54726)],
['C',(0.64866, 1.22109, -1.54261)],
['C',(0.74501, 1.17043, 1.54726)],
['H',(-1.14940, 2.19096, 1.50474)],
['H',(-1.31522, -2.10139, 1.51415)],
['H',(1.14940, -2.19096, 1.50474)],
['H',(3.45997, 0.62992, 1.24374)],
['H',(3.22780, -1.10941, 0.99232)],
['H',(-3.22780, 1.10941, 0.99232)],
['H',(-3.45997, -0.62992, 1.24374)],
['H',(3.22780, 1.10941, -0.99232)],
['H',(3.45997, -0.62992, -1.24374)],
['H',(-3.45997, 0.62992, -1.24374)],
['H',(-3.22780, -1.10941, -0.99232)],
['H',(1.31522, -2.10139, -1.51415)],
['H',(-1.14940, -2.19096, -1.50474)],
['H',(-1.31522, 2.10139, -1.51415)],
['H',(1.14940, 2.19096, -1.50474)],
['H',(1.31522, 2.10139, 1.51415)]]

mol = gto.M(atom = atoms, basis='cc-pvdz') # read the geometry, in A
m = scf.RHF(mol)
m.kernel()

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

fragments,fragment_spins = [[0,16,1],[2,17,3,18],[4,15,31],[5,19,20],[7,23,24],\
                          [6,21,22],[8,25,26],[12,13,29],\
                          [9,14,30],[11,28,10,27]],[0,0,0,0,0,0,0,0,0,0]

thresh   = 1.0e-7
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'cc'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)

