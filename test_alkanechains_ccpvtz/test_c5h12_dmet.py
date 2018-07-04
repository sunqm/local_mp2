import numpy,os
from pyscf import gto,scf,cc,mp
from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb

" C5H12 """
atoms=[
['C',(0.0000000000, 1.2758596583, -0.5229675814)],
['C',(0.0000000000, 0.0000000000, 0.3154943521)],
['C',(0.0000000000, -1.2758596583, -0.5229675814)],
['C',(0.0000000000, -2.5450663093, 0.3242547685)],
['C',(0.0000000000, 2.5450663093, 0.3242547685)],
['H',(0.0000000000, -3.4382954268, -0.3017258094)],
['H',(-0.8771884600, 1.2742283997, -1.1783559671)],
['H',(0.8771884600, 1.2742283997, -1.1783559671)],
['H',(-0.8819324611, -2.5877051429, 0.9666903105)],
['H',(0.8819324611, -2.5877051429, 0.9666903105)],
['H',(0.8776631440, 0.0000000000, 0.9725147596)],
['H',(-0.8776631440, 0.0000000000, 0.9725147596)],
['H',(0.8771884600, -1.2742283997, -1.1783559671)],
['H',(-0.8771884600, -1.2742283997, -1.1783559671)],
['H',(0.0000000000, 3.4382954268, -0.3017258094)],
['H',(-0.8819324611, 2.5877051429, 0.9666903105)],
['H',(0.8819324611, 2.5877051429, 0.9666903105 )]]

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

fragments = [[0,5,6,7],[1,8,9],[2,10,11],[3,12,13],[4,14,15,16]]
fragment_spins = [-1,0,0,0,1]
thresh   = 1.0e-5
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'mp2'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)
