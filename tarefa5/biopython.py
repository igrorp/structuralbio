from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
from Bio.PDB.vectors import calc_angle, calc_dihedral
import seaborn as sns

parser = PDBParser()

structure = parser.get_structure('protease', '6lu7.pdb')



# 1) Calculating bond distances by position

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

# plt.savefig('distances.png', dpi=800)

# 2) Calculating the torsion angles

values = {}

n = 0

for residue in structure.get_residues():

    hetflag, *_ = residue.get_id()

    if hetflag == ' ' and residue.get_resname() != 'GLY':

        nitro = residue['N'].get_vector()
        carbo = residue['C'].get_vector()
        calfa = residue['CA'].get_vector()
        cbeta = residue['CB'].get_vector()
        # angle1 = calc_angle(nitro, carbo, calfa) / 3.14 * 180
        angle1 = calc_angle(nitro, calfa, carbo) / 3.14 * 180 
        angle2 = calc_angle(cbeta, calfa, carbo) / 3.14 * 180
        angle3 = calc_angle(nitro, calfa, cbeta) / 3.14 * 180

        values[n] = [residue.get_resname(), angle1, angle2, angle3]

        n+=1


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

for n, (res, dis1, dis2, dis3) in values.items():

    ax1.scatter(n, dis1, c='red', edgecolors='black')
    ax2.scatter(n, dis2, c='blue', edgecolors='black')
    ax3.scatter(n, dis3, c='green', edgecolors='black')


# plt.savefig('torsionagles.png', dpi=800)


# 3) Calculating diheadral angles by position

residues = []


for residue in structure.get_residues():

    hetflag, *_ = residue.get_id()

    if hetflag == ' ':

        residues.append(residue)


plt.clf()

# data = {'phi':[], 'psi':[]}

for i in range(1, len(residues) - 1):

    nitro = residues[i]['N'].get_vector()
    calfa = residues[i]['CA'].get_vector()
    carbo = residues[i]['C'].get_vector()
    carbomenos1 = residues[i - 1]['C'].get_vector()
    calfamenos1 = residues[i - 1]['CA'].get_vector()
    nitromais1 = residues[i + 1]['N'].get_vector()
    calfamais1 = residues[i + 1]['CA'].get_vector()

    phi = calc_dihedral(calfamenos1, carbomenos1, nitro, calfa)

    phi = calc_dihedral(carbomenos1, nitro, calfa, carbo)

    phi = phi * 180 / 3.14

    psi = calc_dihedral(calfa, carbo, nitromais1, calfamais1)

    psi = calc_dihedral(nitro, calfa, carbo, nitromais1)

    psi = psi * 180 / 3.14

    print(phi, psi)

    # data['phi'].append(phi)
    # data['psi'].append(psi)

    plt.scatter(phi, psi, edgecolors='black')


# sns.kdeplot(data['phi'], data['psi'], cmap='Reds', shade=True, shade_lowest=True, bw=0.1)


plt.savefig('ramachandran.png', dpi=800)