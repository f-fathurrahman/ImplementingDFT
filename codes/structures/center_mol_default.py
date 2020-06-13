#!/home/efefer/miniconda3/bin/python
import sys
import ase.io
from ase.units import Bohr

assert(len(sys.argv) == 2)

atoms = ase.io.read(sys.argv[1])
atoms.set_pbc([True,True,True])
A = 16.0*Bohr
atoms.set_cell([A,A,A])
atoms.center()
atoms.write(sys.argv[1])