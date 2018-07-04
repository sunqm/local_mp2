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
fname="/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen/test_riplinger_copy/DLPNO_CCSD_T-riplinger/210_coord.xyz" #in files

""" 210_coord """
atoms = read_geometry(fname)
mol = gto.M(atom = atoms, basis='cc-pvdz') # read the geometry, in A
m = scf.RHF(mol)
m.kernel()

mycc = cc.CCSD(m)
mycc.kernel()
del mol

bs     = 'dz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs]}
charge = 0
spin   = 0

fragments,fragment_spins = [[0,1,2,3],[4,12,13,14],\
                          [5,15,16,17],[6,18,19,20],\
                          [7,21,22,23],\
                          [8,24,25,26],[9,27,28,29],[10,30,31,32],[11,33,34,35]],[0,-1,1,-1,1,-1,1,-1,1]

thresh   = 1.0e-5
#method = 'cc' uncomment this line and comment the line below to switch to CC
method   = 'cc'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)

