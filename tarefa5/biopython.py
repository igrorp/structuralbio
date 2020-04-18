from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
from Bio.PDB.vectors import calc_angle, calc_dihedral
import seaborn as sns

parser = PDBParser()

structure = parser.get_structure('protease', '6lu7.pdb')

residues = []


for residue in structure.get_residues():

    hetflag, *_ = residue.get_id()

    if hetflag == ' ':

        residues.append(residue)


# 1) Calculating bond distances by position

values = {}

n = 0

for i in range(len(residues) - 1):

    dis1 = residues[i]['CA'] - residues[i]['C']

    dis2 = residues[i]['CA'] - residues[i]['N']

    dis3 = residues[i]['C'] - residues[i+1]['N']

    dis3 = 1.33 if dis3 > 3.0 else dis3

    values[n] = [residues[i].get_resname(), dis1, dis2, dis3]

    print(residues[i].get_resname(), dis3)

    n+=1

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

ax1.hlines(1.52, 0, len(residues))
ax2.hlines(1.45, 0, len(residues))
ax3.hlines(1.33, 0, len(residues))

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_position(('outward', 10))
ax1.set_xticks([])
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_position(('outward', 10))
ax1.set_xticks([])
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_position(('outward', 10))
ax3.spines['left'].set_position(('outward', 10))
ax3.set_xticks([0, 100, 200, 300])

#ax1.set_xticks([], [])

for n, (res, dis1, dis2, dis3) in values.items():

    ax1.scatter(n, dis1, edgecolors='black', c='red')
    ax2.scatter(n, dis2, edgecolors='black', c='blue')
    ax3.scatter(n, dis3, edgecolors='black', c='green')

plt.savefig('distances.png', dpi=800)

# 2) Calculating the torsion angles

plt.clf()

values = {}

n = 0

for residue in residues:

    if residue.get_resname() != 'GLY':

        nitro = residue['N'].get_vector()
        carbo = residue['C'].get_vector()
        calfa = residue['CA'].get_vector()
        cbeta = residue['CB'].get_vector()
        # angle1 = calc_angle(nitro, carbo, calfa) / 3.14 * 180
        angle1 = calc_angle(nitro, calfa, carbo) / 3.14 * 180 
        angle2 = calc_angle(cbeta, calfa, carbo) / 3.14 * 180
        angle3 = calc_angle(nitro, calfa, cbeta) / 3.14 * 180

        if residue.get_resname() == 'SER':

            print(angle1, angle2, angle3)

        values[n] = [residue.get_resname(), angle1, angle2, angle3]

        n+=1


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

for n, (res, dis1, dis2, dis3) in values.items():

    ax1.scatter(n, dis1, c='red', edgecolors='black')
    ax2.scatter(n, dis2, c='blue', edgecolors='black')
    ax3.scatter(n, dis3, c='green', edgecolors='black')

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_position(('outward', 10))
ax1.set_xticks([])
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_position(('outward', 10))
ax1.set_xticks([])
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['left'].set_position(('outward', 10))
ax3.set_xticks([0, 100, 200, 300])

plt.savefig('torsionagles.png', dpi=800)

# Torsion angles by aa type

values = {}

n = 0

for residue in residues:

    if residue.get_resname() != 'GLY':

        nitro = residue['N'].get_vector()
        carbo = residue['C'].get_vector()
        calfa = residue['CA'].get_vector()
        cbeta = residue['CB'].get_vector()
        # angle1 = calc_angle(nitro, carbo, calfa) / 3.14 * 180
        angle1 = calc_angle(nitro, calfa, carbo) / 3.14 * 180 
        angle2 = calc_angle(cbeta, calfa, carbo) / 3.14 * 180
        angle3 = calc_angle(nitro, calfa, cbeta) / 3.14 * 180

        if residue.get_resname() in values:
            
            values[residue.get_resname()].append([angle1, angle2, angle3])

        else:

            values[residue.get_resname()] = [[angle1, angle2, angle3]]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

colors = {
'ALA':'green',
'ARG':'purple',
'ASN':'blue',
'ASP':'red',
'CYS':'yellow',
'GLU':'red',
'GLN':'blue',
'GLY':'green',
'HIS':'purple',
'ILE':'green',
'LEU':'green',
'LYS':'purple',
'MET':'yellow',
'PHE':'green',
'PRO':'green',
'SER':'blue',
'THR':'blue',
'TRP':'green',
'TYR':'green',
'VAL':'green'
}


for aa, data in values.items():

    for (dis1, dis2, dis3) in data:

        ax1.scatter(aa, dis1, c=colors[aa], edgecolors='black')
        ax2.scatter(aa, dis2, c=colors[aa], edgecolors='black')
        ax3.scatter(aa, dis3, c=colors[aa], edgecolors='black')

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
plt.xticks(rotation=45)

plt.savefig('aatorsionagles.png', dpi=800)


# 3) Calculating diheadral angles by position


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

    # data['phi'].append(phi)
    # data['psi'].append(psi)

    plt.scatter(phi, psi, color=colors[residues[i].get_resname()], edgecolors='black')

plt.legend()
plt.savefig('ramachandran.png', dpi=800)


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

    if residue.get_resname() in values:
        
        values[residue.get_resname()].append([phi, psi])

    else:

        values[residue.get_resname()] = [[phi, psi]]


# sns.kdeplot(data['phi'], data['psi'], cmap='Reds', shade=True, shade_lowest=True, bw=0.1)


