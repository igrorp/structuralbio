from Bio.PDB import *

parser = PDBParser()

structure = parser.get_structure('ribozyme', '2n3q.pdb')

sup = Superimposer()

ref = list(structure[0].get_atoms())

for altmodel in structure:

    if altmodel.id != structure[0].id:

        sup.set_atoms(ref, list(altmodel.get_atoms()))

        sup.apply(list(altmodel.get_atoms()))

        print(sup.rotran)

io = PDBIO()
io.set_structure(structure)
io.save('out.pdb')