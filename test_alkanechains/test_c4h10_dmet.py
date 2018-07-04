import numpy
import pyscf
from   pyscf import cc, gto, symm, scf, ao2mo, mp
from   pyscf.tools import localizer
import os
import sys
sys.path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen/')
import orbital_selection_fc as orb

""" C4H10 """
atoms=[
['C',(0.6991898018,  1.8213888862, 0.0000000000)],
['C',(0.7035097406,  0.2955905259, 0.0000000000)],
['C',(-0.7035097406, -0.2955905259, 0.0000000000)],
['C',(-0.6991898018, -1.8213888862, 0.0000000000)],
['H',(1.7121780223,  2.2264196232, 0.0000000000)],
['H',(-1.7121780223, -2.2264196232, 0.0000000000)],
['H',(0.1831112848,  2.2080826508, 0.8815761058)],
['H',(0.1831112848,  2.2080826508, -0.8815761058)],
['H',(-0.1831112848, -2.2080826508, 0.8815761058)],
['H',(-0.1831112848, -2.2080826508, -0.8815761058)],
['H',(1.2467474553, -0.0728437222, -0.8761798942)],
['H',(1.2467474553, -0.0728437222, 0.8761798942)],
['H',(-1.2467474553, 0.0728437222, -0.8761798942)],
['H',(-1.2467474553, 0.0728437222, 0.8761798942 )]]

mol = gto.M(atom=atoms,basis='cc-pvdz')
m   = scf.RHF(mol)
m.kernel()
mm  = cc.CCSD(m)# uncomment this line and comment the line below to switch to CC
#mm = mp.MP2(m)
#mm.kernel()

del mol,m,mm

bs     = 'dz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments = [[0,4,5,6,3,11,12,13],[1,7,8],[2,9,10]]
fragment_spins = [0,0,0]
thresh   = 1.0e-8
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'cc'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)
