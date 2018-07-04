from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb
import numpy,os
from pyscf import gto,scf,cc,mp

""" C2H6 """
atoms=[
['C',( 0.000000, 0.000000, 0.763258)],
['C',(0.000000, 0.000000, -0.763258)],
['H',( 0.000000, 1.016928, 1.157843)],
['H',( -0.880686, -0.508464, 1.157843)],
['H',(  0.880686, -0.508464, 1.157843)],
['H',(  0.000000, -1.016928, -1.157843)],
['H',( -0.880686, 0.508464, -1.157843)],
['H',(  0.880686, 0.508464, -1.157843 )]]

mol = gto.M(atom=atoms,basis='cc-pv5z')
m   = scf.RHF(mol)
m.kernel()
#mm  = cc.CCSD(m)# uncomment this line and comment the line below to switch to CC
mm = mp.MP2(m)
mm.kernel()

del mol,m,mm

bs     = '5z'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments = [[0,2,3,4],[1,5,6,7]]
fragment_spins = [1,-1]
thresh   = 1.0e-8
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'mp2'
nfreeze  = 0
parallel = True

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)
