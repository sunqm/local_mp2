import scipy
import numpy,os
from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb

from pyscf import gto,scf,cc,mp, symm, ao2mo

def read_geometry(fname):
    inf = open(fname,'r')
    x   = inf.readlines()
    x   = x[2:]
    atm = []
    for z in x:
        y = z.split()
        rx,ry,rz = float(y[1]),float(y[2]),float(y[3])
        atm.append([y[0],numpy.asarray([rx,ry,rz])])
    return atm

#------------------------------------------
fname="/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_serial/test_riplinger_copy/DLPNO_CCSD_T-riplinger/anthracene.xyz" #in files

""" anthracene """
atoms = read_geometry(fname)
mol = gto.M(atom = atoms, basis='cc-pvdz') # read the geometry, in A
#m = scf.RHF(mol)
#m.kernel()

# mycc = cc.CCSD(m)
# mycc.kernel()
del mol

bs     = 'dz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments,fragment_spins = [[0,1,2,3,4,5,28,29,30,31],\
                          [8,9,10,11,12,13,34,35,36,37],[14,15,16,17,18,19,38,39,40,41],\
                          [22,23,24,25,26,27,44,45,46,47],\
                          [6,7,20,21,32,33,42,43]],[0,0,0,0,0]

thresh   = 1.0e-7
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'mp2'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)

