import scipy
import numpy,os
from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb

from pyscf import gto,scf,cc,mp, symm, ao2mo, mp

ra = 1.0

atoms=[['N',( 0.000*ra, 1.111*ra,0.000)],
       ['C',(-1.090*ra, 0.282*ra,0.000)],
       ['C',( 1.125*ra, 0.304*ra,0.000)],
       ['N',(-0.750*ra,-0.994*ra,0.000)],
       ['C',( 0.638*ra,-0.988*ra,0.000)],
       ['H',(-0.015*ra, 2.123*ra,0.000)],
       ['H',(-2.111*ra, 0.668*ra,0.000)],
       ['H',( 2.135*ra, 0.712*ra,0.000)],
       ['H',( 1.210*ra,-1.916*ra,0.000)]]

mol = gto.M(atom = atoms, basis='cc-pvdz') # read the geometry, in A
m = scf.RHF(mol)
m.kernel()

#mycc = cc.CCSD(m)
#mycc.kernel()

#mymp2 = mp.MP2(m)
#mymp2.kernel()
del mol

bs     = 'dz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs, 'N': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs], 'N': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments = [[0,4,8],[3,7],[2,6],[1,5]]
fragment_spins = [0,0,0,0]

thresh   = 1.0e-6
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'mp2'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)
