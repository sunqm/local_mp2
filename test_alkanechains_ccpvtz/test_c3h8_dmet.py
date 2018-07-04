import numpy
import pyscf
from   pyscf import cc, gto, symm, scf, ao2mo, mp
from   pyscf.tools import localizer
import os
import sys
sys.path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen/')
import orbital_selection_fc as orb

""" C3H8 """
atoms=[
['C',( 0.000000, 1.265790, -0.260971)],
['H',( -0.881481, 1.298774, -0.905267)],
['H',( 0.881481, 1.298774, -0.905267)],
['H',( 0.000000, 2.166499, 0.354619)],
['C',( 0.000000, 0.000000, 0.591729)],
['H',( 0.875492, 0.000000, 1.246553)],
['H',( -0.875492, 0.000000, 1.246553)],
['C',( 0.000000 ,-1.265790, -0.260971)],
['H',( 0.881481, -1.298774, -0.905267)],
['H',( 0.000000 ,-2.166499 ,0.354619)],
['H',( -0.881481, -1.298774 ,-0.905267 )]]

mol = gto.M(atom=atoms,basis='cc-pvtz')
m   = scf.RHF(mol)
m.kernel()
#mm  = cc.CCSD(m)# uncomment this line and comment the line below to switch to CC
mm = mp.MP2(m)
mm.kernel()

del mol,m,mm

bs     = 'tz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments = [[0,1,2,3],[4,5,6],[7,8,9,10]]
fragment_spins = [-1,0,1]
thresh   = 1.0e-8
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'mp2'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)
