from Bio import PDB
from numpy import array
from math import pi

parser = PDB.PDBParser()            # Estrutura para parsear arquivo PDB
io = PDB.PDBIO()                    # Salva em variavel
struct = parser.get_structure('253l',"253l.pdb")        # Pega estrutura escolhida
                                                        # Aqui depende do seu pdb


arq = open("movie_nmr.pdb","w")         # Abre arquivo nmr de output, que no pymol passa como movie
                                        # Com isso, só é necessário abrir o movie_nmr.pdb no pymol

num_total = 0                           # Numero total de atomos
coord_x = 0                             # Soma total de coordenadas x
coord_y = 0                             # Soma total de coordenadas y
coord_z = 0                             # Soma total de coordenadas z

for residue in struct.get_residues():
    for atom in residue:
        num_total += 1
        coord_x += atom.get_coord()[0]
        coord_y += atom.get_coord()[1]
        coord_z += atom.get_coord()[2]   # Aqui pega as coordenadas de cada atomo em cada eixo e faz soma total

coord_x = coord_x/num_total
coord_y = coord_y/num_total
coord_z = coord_z/num_total              # Torna o somatório na média em cada eixo, ou seja, centróide


z = 0.1                                     # Variável que salva quantos radianos será girado o sistema
                                            # Quanto menor, mais suave será
for i in range(0,50):
    translation = array((0, 0 , 0), 'f')    # Matriz de translação. Se quiser transladar, mudar valores do array

    if i != 0:                              # Se nao eh a primeira estrutura
        struct = parser.get_structure('teste{}'.format(i-1), "teste{}.pdb".format(i-1))  # Pega o frame anterior
        for residue in struct.get_residues():                   
            for atom in residue:
                coord_x += atom.get_coord()[0]
                coord_y += atom.get_coord()[1]
                coord_z += atom.get_coord()[2]

        coord_x = coord_x / num_total
        coord_y = coord_y / num_total
        coord_z = coord_z / num_total                                                   # E faz o mesmo que foi feito antes
                                                                                        # Gerando novo centroide

    for residue in struct.get_residues():           # Para todo residuo na estrutura
        #res_name = residue.get_resname()            # Pega nome do residuo
        #if res_name == "BGC":                       # Se residuo eh glicose
        for atom in residue:                    # Para todos os atomos

            #rotation_matrix = PDB.rotmat(PDB.Vector(0,0,z), PDB.Vector(0,0,0))  # Matriz de rotacao
            rotation_matrix = PDB.rotmat(PDB.Vector(z,z,z), PDB.Vector(coord_x,coord_y,coord_z))               # Matriz de rotacao em si mesmo
            #rotation_matrix = PDB.rotmat(PDB.Vector(1, 1, 1), PDB.Vector(0.5,0,0))               # Matriz para rotar na central
            atom.transform(rotation_matrix,translation)                                           # Aplica transformacao

    io = PDB.PDBIO()
    io.set_structure(struct)
    io.save("teste{}.pdb".format(i))                               # Salva estrutura do frame em .pdb


    # Rotina para criar arquivo NMR que pode ser aberto no pymol ja como um movie

    with open("teste{}.pdb".format(i),"r") as file:
        arq.write("MODEL\t {}\n".format(i))
        for lines in file:
            if not "END" in lines:
                arq.write(lines)
        arq.write("ENDMDL")
