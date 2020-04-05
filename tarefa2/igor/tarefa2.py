
file = '6ofy.pdb'


class Atom():

    def __init__(self, classifier, natom, atomname, altlocation, resname, chain, nresidue, resinsertion, x, y, z, occ, tfactor, elementsymbol, atomcharge):

        self.classifier = classifier
        self.natom = int(natom)
        self.atomname = atomname
        self.altlocation = altlocation
        self.resname = resname
        self.chain = chain
        self.nresidue = int(nresidue)
        self.resinsertion = resinsertion
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.occ = float(occ)
        self.tfactor = float(tfactor) 
        self.elementsymbol = elementsymbol
        self.atomcharge = atomcharge


class PDB():

    def __init__(self, filename):

        self.atoms = []

    # Instanciating only the ATOM lines for now

        with open(filename) as file:

            for number, line in enumerate(file.readlines()):

                classifier = line[0:6]

                if classifier == 'ATOM  ' or classifier == 'HETATM':
                    
                    try:

                        classifier = line[0:6]
                        natom = line[6:11]
                        atomname = line[12:16]
                        altlocation = line[16:17]
                        resname = line[17:20]
                        chain = line[21:22]
                        nresidue = line[22:26]
                        resinsertion = line[26:27]
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        occ = line[54:60]
                        tfactor = line[60:66]
                        elementsymbol = line[76:78]
                        atomcharge = line[78:80]
                                           
                        self.atoms.append(Atom(classifier, natom, atomname, altlocation, resname, chain, nresidue, resinsertion, x, y, z, occ, tfactor, elementsymbol, atomcharge))

                    except:

                        print(f'Could not parse line {number}')

    
    def dump_to_file(self, filename):

        '''Saves the current structure info into a .pdb file'''

        with open(filename, 'w') as newfile:

            for atom in self.atoms:

                newfile.write(
                    "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(
                        atom.classifier,
                        atom.natom,
                        atom.atomname,
                        atom.altlocation,
                        atom.resname,
                        atom.chain,
                        atom.nresidue,
                        atom.resinsertion,
                        atom.x,
                        atom.y,
                        atom.z,
                        atom.occ,
                        atom.tfactor,
                        atom.elementsymbol,
                        atom.atomcharge
                    )
                )

    def move_atoms(self, dx, dy, dz, resname=None):

        for atom in self.atoms:
            
            if not resname or atom.resname == resname:

                atom.x+=float(dx)
                atom.y+=float(dy)
                atom.z+=float(dz)


    def create_model(self, n, dx, dy, dz, resname=None, filename='model.pdb'):

        with open(filename, 'w') as newfile:
            
            for i in range(n):

                newfile.write(f'MODEL        {i+1}\n')
            
                for atom in self.atoms:

                    if not resname or resname == atom.resname:

                        newfile.write(
                        "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(
                            atom.classifier,
                            atom.natom,
                            atom.atomname,
                            atom.altlocation,
                            atom.resname,
                            atom.chain,
                            atom.nresidue,
                            atom.resinsertion,
                            atom.x + (i * dx),
                            atom.y + (i * dy),
                            atom.z + (i * dz),
                            atom.occ,
                            atom.tfactor,
                            atom.elementsymbol,
                            atom.atomcharge
                        ))

                newfile.write('ENDMDL\n')


file = '6ofy.pdb'
prostaglandin = PDB(file)

prostaglandin.create_model(50, -0.5, 0, 0, filename='modelx.pdb', resname='COH')
prostaglandin.create_model(50, 0, -0.5, 0, filename='modely.pdb', resname='COH')
prostaglandin.create_model(50, 0, 0, -0.5, filename='modelz.pdb', resname='COH')

# prostaglandin.move_atoms(-5, -5, -5, 'COH')
# prostaglandin.dump_to_file('test.pdb')

