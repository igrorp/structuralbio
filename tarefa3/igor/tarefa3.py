import numpy as np
import math

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
        self.coords = np.array([[float(x), float(y), float(z)]]).transpose()
        #print(self.coords)
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
                        atom.coords[0][0],
                        atom.coords[1][0],
                        atom.coords[2][0],
                        atom.occ,
                        atom.tfactor,
                        atom.elementsymbol,
                        atom.atomcharge
                    )
                )


    def move_atoms(self, dx, dy, dz, resname=None):

        for atom in self.atoms:
            
            if not resname or atom.chain == resname:

                atom.coords[0][0] = float(dx) - atom.coords[0][0]
                atom.coords[1][0] = float(dy) - atom.coords[1][0]
                atom.coords[2][0] = float(dz) - atom.coords[2][0]


    def translate_model(self, n, dx, dy, dz, resname=None, filename='model.pdb'):

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
                            atom.coords[0][0] + (i * dx),
                            atom.coords[1][0] + (i * dy),
                            atom.coords[2][0] + (i * dz),
                            atom.occ,
                            atom.tfactor,
                            atom.elementsymbol,
                            atom.atomcharge
                        ))

                newfile.write('ENDMDL\n')

    
    def centroid(self, chain=None):

        x, y, z, len = 0, 0, 0, 0

        for atom in self.atoms:

            if not chain or chain == atom.chain:

                x+=atom.coords[0][0]
                y+=atom.coords[1][0]
                z+=atom.coords[2][0]

                len+=1

        print(x / len, y / len, z / len)

        self.move_atoms(x / len, y / len, z / len, resname=chain)


    def rotate(self, dx=0, dy=0, dz=0, steps=1, chain=None, filename='model.pdb'):

        cosx = math.cos(dx * math.pi / 180 / steps)
        cosy = math.cos(dy * math.pi / 180 / steps)
        cosz = math.cos(dz * math.pi / 180 / steps)
        sinx = math.sin(dx * math.pi / 180 / steps)
        siny = math.sin(dy * math.pi / 180 / steps)
        sinz = math.sin(dz * math.pi / 180 / steps)

        rotmatrix = np.array([
            #[cosz * cosy, cosz * siny * sinx - sinz * cosx, cosz * siny * sinx + sinz * cosx],
            [cosz * cosy, cosz * siny * sinx - sinz * cosx, cosz * siny * cosx + sinz * sinx],
            #[sinz * cosy, sinz * siny * sinx + cosz * cosx, sinz * siny * sinx + cosz * cosx],
            [sinz * cosy, sinz * siny * sinx + cosz * cosx, sinz * siny * cosx - cosz * sinx],
            [-1 * siny, cosy * sinx, cosy * cosx]
        ])

        with open(filename, 'w') as newfile:

            for n in range(steps):

                newfile.write(f'MODEL        {n+1}\n')

                for atom in self.atoms:

                    if not chain or chain == atom.chain:

                        atom.coords = rotmatrix.dot(atom.coords)

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
                            atom.coords[0][0],
                            atom.coords[1][0],
                            atom.coords[2][0],
                            atom.occ,
                            atom.tfactor,
                            atom.elementsymbol,
                            atom.atomcharge
                        ))


                newfile.write('ENDMDL\n')



trna = PDB('1g59.pdb')
trna.centroid(chain='B')
#trna.rotate(360, 0, 0, 50, chain='B', filename='modelX.pdb')
#trna.rotate(0, 360, 0, 50, chain='B', filename='modelY.pdb')
trna.rotate(0, 0, 360, 50, chain='B', filename='modelZ.pdb')
