import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
import numpy as np
import argparse

# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type = str, default='martini_graphene', help = 'Name of the output, default = martini_graphene')

args = parser.parse_args()

filename = args.filename


while True:
    rows = int(input("Enter the number of rows, make sure that it is odd, and greater than 3: "))
    if (rows >=3) and (rows % 2 !=0):
    # the correct input was entered, so we break out of the loop
        break
    else:
    # the incorrect input was entered, so we continue the loop
        continue


while True:
    columns = int(input("Enter the number of beads across a row, make sure that it is a multiple of 3: "))
    if (columns % 3 == 0):
        columns = columns - 1
        break
    else:
    # the incorrect input was entered, so we continue the loop
        continue



#----------------#
# Structure File #
#----------------#



positions = []
dist = 0 #setting y = 0 in the beginning



for i in range(rows):
    if i % 2 == 0: # if the row is indexed even number
        
        if (i == rows-1) or (i == 0): # special treatment for the first and the last row
            x = np.linspace(0, (columns-1) * 2.56, num = columns) #AP-> tn = a+(n-1)d
            d_array = np.array(np.arange(3, columns+1, 3)-1) # array for deleting the elements from the first and the last row which are virtual sites
            x = np.delete(x, d_array) #deleting those virtual sites from the first and the last row
        else: # if the row is even. but are not the first and the last one
            x = np.linspace(0, (columns-1) * 2.56, num = columns) #AP-> tn = a+(n-1)d
            
            
        
        for k in range(len(x)): # Now looping every x-coordinate
            positions.append([x[k], dist, 0]) # Adding the x,y,z coordinate for the position
            
            
        
        dist+=2.217 #Going 2.17 Angstroms down/up along y

        
    else: # if the row is indexed odd number
        x = np.linspace(-1.28, -1.28+columns*2.56, num = columns+1) #AP-> tn = a+(n-1)d

        for k in range(len(x)):
            positions.append([x[k], dist, 0]) 
        
        
        dist+=2.217 #Going 2.17 Angstroms down/up along y
        

        

    

w = mda.Universe.empty(n_atoms = len(positions), trajectory = True) # Creating an empty Universe with atoms, and trajectory is True for writing positions
w.atoms.positions = positions
w.atoms.write('input.gro')


u = mda.Universe(filename+'.gro')
c = np.unique(u.atoms.positions[:,1]) # finding the unique y-coordinate so that we can loop through each row of atoms
u.atoms.masses = 36    # Since it's a TC5 bead, we use a mass of 36 for each


'''
Identifying the indices of the virtual sites, and setting their mass to zero
'''
def virtual_site(universe):
    ''' Return the serial or (1-based) index for virtual sites, and also sets the mass for the virtual-site to 0'''
    for j in range(len(c)): # Looping through each row, defined by the y-coordinate
        if (j !=0) and (j !=len(c)-1):
            if j % 2 !=0:
                group = u.atoms[u.atoms.positions[:,1] == c[j]]
                gr = np.arange(1, len(group), 3)
                print(gr)
                for k in gr:
                    u.atoms[group[k].index].mass = 0
    
    return u.atoms[u.atoms.masses == 0].indices + 1

    
for j in range(len(c)): # Looping through each row, defined by the y-coordinate
    if (j !=0) and (j !=len(c)-1):
        if j % 2 !=0:
            group = u.atoms[u.atoms.positions[:,1] == c[j]]
            gr = np.arange(1, len(group), 3)
            print(gr)
            for k in gr:
                u.atoms[group[k].index].mass = 0
            
        else:
            group = u.atoms[u.atoms.positions[:,1] == c[j]]
            gr = np.arange(2, len(group), 3)
            for k in gr:
                u.atoms[group[k].index].mass = 0
                
                


def hexagon(universe):
    
    ''' Identifying the vertics of hexagon around a virtual site '''
    
    b = u.atoms[u.atoms.masses == 0].indices
    hexagon_indices = []
    for i in b:
        empty = []
        for j in u.atoms.indices:
            if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                empty.append(j)
        hexagon_indices.append(empty)
    return hexagon_indices



def bonds(universe):
    ''' Returns the list of pair of indices for which bonds are defined '''
    list_of_bonds = []
    for element in hexagon(u):
        for i in element:
            for j in element:
                if i < j:
                    if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                        list_of_bonds.append((i,j))
    return list_of_bonds
                        
                        
def angles(universe):
    ''' Returns the list of triplet of indices for which an angle is defined '''
    list_of_angles = []
    for element in hexagon(u):
        for i in range(len(element)):
            list_of_angles.append((element[i-1], element[i], element[(i+1) % len(element)])) # % operator is used to wrap the indices around the end of the list, so that the triplets at the beginning and end of the list include the first and last elements
    
    
            
            
    return list_of_angles



def virtual_sites(universe):
    
    ''' Identifying the vertics of hexagon around a virtual site '''
    
    b = u.atoms[u.atoms.masses == 0].indices
    hexagon_indices = []
    for i in b:
        empty = []
        
        for j in u.atoms.indices:
            if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                empty.append(u.atoms[j])
        empty = sorted(empty, key = lambda x: x.position[1])
        empty = [i.index for i in empty]
        data = empty[:2] + empty[-2:]
        data.append(u.atoms[i].index)
        data = data[::-1]
        hexagon_indices.append(data)
        
    return hexagon_indices



def virtual_sites(universe):
    
    ''' Identifying the vertics of hexagon around a virtual site '''
    
    b = u.atoms[u.atoms.masses == 0].indices
    hexagon_indices = []
    for i in b:
        empty = []
        empty.append(i)
        for j in u.atoms.indices:
            if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                empty.append(j)
        
        hexagon_indices.append(empty)
    return hexagon_indices


#---------------#
# Topology File #
#---------------#


# Open the file for writing

topology_file = open(filename+".itp", 'w')

# Variables



# Header

topology_file.write( "; \n;  Graphene topology\n; for the Martini3 force field\n;\n; created by martini3-graphene-topology.py\n;\n" )
topology_file.write( "; Roshan Shrestha\n; CNRS\n;\n\n" )

topology_file.write("[ moleculetype ]\n")
topology_file.write("; molname	 nrexcl\n")
topology_file.write("  GRA           1")



# Atoms

topology_file.write( "\n[ atoms ]\n" )
topology_file.write( "; nr	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n" )

#for i in range(1, numatoms+1):



# Bonds

topology_file.write( "\n[ bonds ]\n" )
topology_file.write( "; i	 j	  funct	 length	 kb\n" )










# Angles

topology_file.write( "\n[ angles ]\n" )
topology_file.write( "; i	 j	 k	 funct	 angle	 force_k\n" )





# Improper Dihedrals

topology_file.write( "\n[ dihedrals ]\n" )
topology_file.write( "; i	 j	 k	 l     funct	 ref.angle     force_k\n" )




# Virtual sites

topology_file.write( "\n[ virtual_sitesn ]\n" )
topology_file.write( "; site	 funct	 constructing atom indices\n" )
for i in virtual_sites(u):
    topology_file.write(f"  {i[0]+1:<3}     1     {i[1]+1:<3}     {i[2]+1:<3}     {i[3]+1:<3}    {i[4]+1:<3}\n")

# Exclusions
topology_file.write( "\n[ exclusions ]\n" )
for i in exclusions(u):
    topology_file.write(f"  {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     {i[3]+1:<3}    {i[4]+1:<3}     {i[5]+1:<3}     {i[6]+1:<3}\n")

topology_file.close()




                        
        
        

        


    
        
                

            
        



        


        
