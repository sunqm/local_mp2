import scipy
import numpy,os
from sys import path
path.append('/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_parallel_ccsdt_frozen')
import orbital_selection_fc as orb
from   pyscf.tools import localizer
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
fname="/home/yuliya/DMET/alkanechains_scf_ccsd/dmet_serial/test_riplinger_copy/DLPNO_CCSD_T-riplinger/penicilline.xyz" #in files

""" penicilline """
atoms = read_geometry(fname)
mol = gto.M(atom = atoms, basis='cc-pvdz') # read the geometry, in A
#m = scf.RHF(mol)
#m.kernel()

# mycc = cc.CCSD(m)
# mycc.kernel()
del mol

bs     = 'dz'
basis  = {'C': 'cc-pv'+bs, 'H': 'cc-pv'+bs, 'O': 'cc-pv'+bs, 'N': 'cc-pv'+bs, 'S': 'cc-pv'+bs,}
shells = {'C': ['sto-6g','cc-pv'+bs], 'H': ['sto-6g','cc-pv'+bs], 'O': ['sto-6g','cc-pv'+bs], 'N': ['sto-6g','cc-pv'+bs], 'S': ['sto-6g','cc-pv'+bs]}

thresh   = 1.0e-5
charge = 0
spin   = 0

fragments,fragment_spins = [[18,19,20,21,22,23,37,38,39,30,41],[15,17,32,36],[16,14,9,27],[0,1,2,3,4,5,6,8,24,25,26],[10,11,28,29,30,31,32,33],[7,12,13,34]],[0,0,0,0,0,0]

method   = 'cc'
nfreeze  = 0
parallel = False

orb.DMET_wrap(atoms,basis,charge,spin,fragments,fragment_spins,shells,nfreeze,method,thresh,parallel)


