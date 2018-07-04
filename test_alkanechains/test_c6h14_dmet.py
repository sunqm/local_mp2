import numpy
import pyscf
from   pyscf import cc, gto, symm, scf, ao2mo, mp
from   pyscf.tools import localizer
import os
import sys
sys.path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen/')
import orbital_selection_fc as orb

""" C6H14 """
atoms=[
['C',(3.1987314968, -0.2104570612, 0.0000000000)],
['H',(3.2763256850, -0.8506634635, 0.8816038652)],
['H',(4.0552248738, 0.4654264605, 0.0000000000)],
['H',(3.2763256850, -0.8506634635, -0.8816038652)],
['C',(1.8773135249, 0.5535425093, 0.0000000000)],
['H',(1.8289548319, 1.2081298850, 0.8760991121)],
['H',(1.8289548319, 1.2081298850, -0.8760991121)],
['C',(0.6649585374, -0.3739765481, 0.0000000000)],
['H',(0.7125760695, -1.0301702097, 0.8769078936)],
['H',(0.7125760695, -1.0301702097, -0.8769078936)],
['C',(-0.6649585374, 0.3739765481, 0.0000000000)],
['H',(-0.7125760695, 1.0301702097, -0.8769078936)],
['H',(-0.7125760695, 1.0301702097, 0.8769078936)],
['C',(-1.8773135249, -0.5535425093, 0.0000000000)],
['H',(-1.8289548319, -1.2081298850 ,0.8760991121)],
['H',(-1.8289548319, -1.2081298850, -0.8760991121)],
['C',(-3.1987314968, 0.2104570612, 0.0000000000)],
['H',(-3.2763256850, 0.8506634635, -0.8816038652)],
['H',(-4.0552248738, -0.4654264605, 0.0000000000)],
['H',(-3.2763256850, 0.8506634635, 0.8816038652)]]
basis ='cc-pvdz'

mol = gto.M(atom=atoms,basis='cc-pvdz')
m   = scf.RHF(mol)
m.kernel()
#mm  = cc.CCSD(m)# uncomment this line and comment the line below to switch to CC
mm = mp.MP2(m)
mm.kernel()

del mol,m,mm

bs     = 'dz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments = [[0,1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14,15],[16,17,18,19]]
fragment_spins = [-1,0,0,0,0,1]
thresh   = 1.0e-8
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'mp2'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)
