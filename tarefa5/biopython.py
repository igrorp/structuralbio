from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
from Bio.PDB.vectors import calc_angle, calc_dihedral

parser = PDBParser()

structure = parser.get_structure('protease', '6lu7.pdb')

# Calculating bond distances

values = {}

n = 0

for residue in structure.get_residues():

    hetflag, *_ = residue.get_id()

    if hetflag == ' ':

        dis1 = residue['CA'] - residue['C']

        dis2 = residue['CA'] - residue['N']

        dis3 = residue['C'] - residue['N']

        values[n] = [residue.get_resname(), dis1, dis2, dis3]

        n+=1

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

for n, (res, dis1, dis2, dis3) in values.items():

    ax1.scatter(n, dis1, edgecolors='black', c='red')
    ax2.scatter(n, dis2, edgecolors='black', c='blue')
    ax3.scatter(n, dis3, edgecolors='black', c='green')

plt.savefig('distances.png', dpi=800)

# Calculating the torsion angles

values = {}

n = 0

for residue in structure.get_residues():

    hetflag, *_ = residue.get_id()

    if hetflag == ' ' and residue.get_resname() != 'GLY':

        nitro = residue['N'].get_vector()
        carbo = residue['C'].get_vector()
        calfa = residue['CA'].get_vector()
        cbeta = residue['CB'].get_vector()
        angle1 = calc_angle(nitro, carbo, calfa) / 3.14 * 180 
        angle2 = calc_angle(carbo, cbeta, calfa) / 3.14 * 180
        angle3 = calc_angle(cbeta, calfa, nitro) / 3.14 * 180

        values[n] = [residue.get_resname(), angle1, angle2, angle3]

        n+=1


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

for n, (res, dis1, dis2, dis3) in values.items():

    ax1.scatter(n, dis1, c='red', edgecolors='black')
    ax2.scatter(n, dis2, c='blue', edgecolors='black')
    ax3.scatter(n, dis3, c='green', edgecolors='black')


plt.savefig('torsionagles.png', dpi=800)

        