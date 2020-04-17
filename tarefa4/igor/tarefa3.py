import numpy as np
import math
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
        self.occ = float(occ)
        self.tfactor = float(tfactor) 
        self.elementsymbol = elementsymbol
        self.atomcharge = atomcharge

    def rotate(self, dx=0, dy=0, dz=0, chain=None):

        cosx = math.cos(dx * math.pi / 180)
        cosy = math.cos(dy * math.pi / 180)
        cosz = math.cos(dz * math.pi / 180)
        sinx = math.sin(dx * math.pi / 180)
        siny = math.sin(dy * math.pi / 180)
        sinz = math.sin(dz * math.pi / 180)

        rotmatrix = np.array([
            [cosz * cosy, cosz * siny * sinx - sinz * cosx, cosz * siny * cosx + sinz * sinx],
            [sinz * cosy, sinz * siny * sinx + cosz * cosx, sinz * siny * cosx - cosz * sinx],
            [-1 * siny, cosy * sinx, cosy * cosx]
        ])

        if not chain or chain == self.chain:

            #self.coords = rotmatrix.dot(self.coords)

            newcoords = rotmatrix.dot(self.coords)

            return Atom(
                self.classifier,
                self.natom,
                self.atomname,
                self.altlocation,
                self.resname,
                self.chain,
                self.nresidue,
                self.resinsertion,
                newcoords[0][0],
                newcoords[1][0],
                newcoords[2][0],
                self.occ,
                self.tfactor,
                self.elementsymbol,
                self.atomcharge
            )

    def translate(self, dx=0, dy=0, dz=0, chain=None):

        transarray = np.array([[dx, dy, dz]]).transpose()

        if not chain or chain == self.chain:

            self.coords = transarray + self.coords



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


    def rotation_model(self, dx=0, dy=0, dz=0, steps=1, chain=None, filename='model.pdb'):

        cosx = math.cos(dx * math.pi / 180 / steps)
        cosy = math.cos(dy * math.pi / 180 / steps)
        cosz = math.cos(dz * math.pi / 180 / steps)
        sinx = math.sin(dx * math.pi / 180 / steps)
        siny = math.sin(dy * math.pi / 180 / steps)
        sinz = math.sin(dz * math.pi / 180 / steps)

        rotmatrix = np.array([
            [cosz * cosy, cosz * siny * sinx - sinz * cosx, cosz * siny * cosx + sinz * sinx],
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



class Model():

    def __init__(self, filename):

        self.models = []

        # A list of PDB objects containing each one of the models inside de file

        with open(filename) as userfile:

            currmodel = []

            for number, line in enumerate(userfile.readlines()):

                classifier = line[0:6]

                if classifier == 'MODEL ':

                    self.models.append(currmodel)

                    currmodel = []

                elif classifier == 'ATOM  ':

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
                                        
                        currmodel.append(
                            Atom(
                                classifier,
                                natom,
                                atomname,
                                altlocation,
                                resname,
                                chain,
                                nresidue,
                                resinsertion,
                                x,
                                y,
                                z,
                                occ,
                                tfactor,
                                elementsymbol,
                                atomcharge)
                        )

                    except:

                        print(f'Could not parse line {number}')
        
        self.models.pop(0)
    
    def rmsd(self, ref, cmp):#ref=0, cmp=1):

        sum, len = 0, 0

        #refmodel = self.models[ref]

        #altmodel = self.models[cmp]

        refmodel = ref
        altmodel = cmp

        # if len(refmodel) != len(altmodel):

        #     raise Exception('Models have different number of atoms')

        for refatom, altatom in zip(refmodel, altmodel):

            atomsum = (refatom.coords - altatom.coords) ** 2

            sum += atomsum.sum()

            len+=1

        #refmodel = self.models[ref]

        return f"{math.sqrt(sum / len):.20f}"
                
# ribozyme = Model('2n3q.pdb')


ribozyme = Model('2n3q.pdb')

print('Original RMSD -->', ribozyme.rmsd(ribozyme.models[0], ribozyme.models[1]))

xref, yref, zref, reflen = 0, 0, 0, 0

for atom in ribozyme.models[0]:

    xref+=atom.coords[0][0]
    yref+=atom.coords[1][0]
    zref+=atom.coords[2][0]

    reflen+=1

xalt, yalt, zalt, altlen = 0, 0, 0, 0

for atom in ribozyme.models[1]:

    xalt+=atom.coords[0][0]
    yalt+=atom.coords[1][0]
    zalt+=atom.coords[2][0]

    altlen+=1

for atom in ribozyme.models[1]:

    atom.translate((xalt - xref) / reflen, (yalt - yref) / reflen, (zalt - zref) / reflen)

print('After centroid translation -->', ribozyme.rmsd(ribozyme.models[0], ribozyme.models[1]))

currval = float(ribozyme.rmsd(ribozyme.models[0], ribozyme.models[1]))
a = currval
newval = 100

counter = 0

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

dict = {'x':[], 'y':[], 'z':[], 'rmsd':[]}

for _ in range(1000):

    x = (random.random() - 0.5) * 10
    y = (random.random() - 0.5) * 10
    z = (random.random() - 0.5) * 10

    newatoms = []

    for atom in ribozyme.models[1]:

        newatoms.append(atom.rotate(x, y, z))
    
    newval = float(ribozyme.rmsd(ribozyme.models[0], newatoms))

    dict['x'].append(x)
    dict['y'].append(y)
    dict['z'].append(z)
    dict['rmsd'].append(a / newval)

    if newval < currval:

        currval = newval

print(f'--> {x}, {y}, {z} --> {currval}')


ax.scatter(dict['x'], dict['y'], dict['z'], c=dict['rmsd'], cmap=plt.get_cmap('Greens'))
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
#plt.colorbar()
plt.savefig('test2.png', dpi=800)



