import numpy
from   pyscf import cc, gto, symm, scf, ao2mo, mp
import scipy
import sys

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
basis='dz'

print (read_geometry(fname))

mol = gto.M(
    atom = read_geometry(fname),                          # read the geometry, in Angstrom, from the given file
    basis = {'C': 'cc-pv'+basis, 'H': 'cc-pv'+basis},     # define the heavy-augmented basis
    verbose = 4,
    max_memory = 24000)
m = scf.RHF(mol)
m.kernel()

#mycc = cc.CCSD(m)
#mycc.kernel()

mymp2 = mp.MP2(m)
mymp2.kernel()
