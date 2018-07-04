from pyscf import gto, scf, lo
from pyscf.tools import localizer

mol = gto.Mole()
mol.atom = '''
     He   0.    0.     0.2
     H    0.   -0.5   -0.4
     H    0.    0.5   -0.4
  '''
mol.basis = 'sto-3g'
mol.build()
mf = scf.RHF(mol).run()

print mf.mo_coeff

loc  = localizer.localizer( mol, mf.mo_coeff, 'boys' )
loc.verbose = 5
print( "loc object")
print (loc.coeff)
new_coeff = loc.optimize(threshold=1e-5)
#loc = lo.Boys(mol, mf.mo_coeff).kernel()
print (new_coeff)
#print("localized = ")
#print( loc)

