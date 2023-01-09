import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
import numpy as np
import argparse
import math

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
w.atoms.write(filename+'.gro')


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
                
                for k in gr:
                    u.atoms[group[k].index].mass = 0
    
    return u.atoms[u.atoms.masses == 0].indices + 1

    
for j in range(len(c)): # Looping through each row, defined by the y-coordinate
    if (j !=0) and (j !=len(c)-1):
        if j % 2 !=0:
            group = u.atoms[u.atoms.positions[:,1] == c[j]]
            gr = np.arange(1, len(group), 3)
            
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
    bds = []
    for element in hexagon(u):
        for i in element:
            for j in element:
                if i < j:
                    if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                        bds.append((i,j))
                        
    for bond in bds:
        if bond not in list_of_bonds:
            list_of_bonds.append(bond)
                        
    
    return list_of_bonds
                        
                        
def angles(universe):
    ''' Returns the list of triplet of indices for which an angle is defined '''
    list_of_angles = []
    dd = []
    for hex in hexagon(u):
        for i in hex:
            for j in hex:
                for k in hex:
                    if (119 <= math.degrees(calc_angles(u.atoms[i].position, u.atoms[j].position, u.atoms[k].position)) <= 121):
                        dd.append([i,j,k])
                        
    for angle in dd:
        if angle[0] > angle[-1]:
            angle[0], angle[-1] = angle[-1], angle[0]
            
    for val in dd:
        if val not in list_of_angles:
            list_of_angles.append(val)
            
    return list_of_angles



def virtual_sites(universe):
    
    ''' Returns the list of virtual site, and its corresponding four beads from which it is constructed '''
    
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

for i in virtual_sites(u):
    for j in i[1:]:
        u.atoms[j].mass +=9



def exclusions(universe):
    
    ''' Returns the list of virtual site and its neighbouring 6 beads for exclusions'''
    
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


hexagons = hexagon(u)
def get_angles(hexagons):
    angles = []
    for hex in hexagons:
        angle1 = [hex[2], hex[4], hex[0], hex[5]]
        angle2 = [hex[4], hex[0], hex[5], hex[1]]
        angle3 = [hex[0], hex[5], hex[1], hex[3]]
        angles.append(angle1)
        angles.append(angle2)
        angles.append(angle3)
    return angles



def get_between_hexagons_angles(hexagons, beads_per_row):
    def get_common_side(hexa1, hexa2):
        common = list(set(hexa1).intersection(set(hexa2)))
        common.sort()
        return common
    
    angles = []
    gap_for_above = int(beads_per_row / 2)
    
    for i in range(len(hexagons)):
        hexa = hexagons[i]
        try:
            hexa_after = hexagons[i+2]
        except IndexError:
            hexa_after = None
        try: 
            hexa_above = hexagons[i+gap_for_above]
        except IndexError:
            hexa_above = None
        
        try: 
            hexa_before = hexagons[i - 1]
        except IndexError:
            hexa_before = None
            
        
        if hexa_above == hexa_before: 
            hexa_before = None
        
        if hexa_after:
            side_after = get_common_side(hexa, hexa_after)
            if side_after:
                idx1 = hexa.index(side_after[0])
                idx2 = hexa_after.index(side_after[1])
                angle = [hexa[idx1 - 2], side_after[0], side_after[1], hexa_after[idx2 + 2]]
                angles.append(angle)
            
        
        if hexa_above: 
            side_above = get_common_side(hexa, hexa_above)
            if side_above: 
                idx1 = hexa.index(side_above[0])
                idx2 = hexa_above.index(side_above[1])
                angle = [hexa[idx1 - 4], side_above[0], side_above[1], hexa_above[idx2 + 4]]
                angles.append(angle)  
                
        if hexa_before:
            side_before = get_common_side(hexa, hexa_before)
            if side_before:
                idx1 = hexa_before.index(side_before[0])
                idx2 = hexa.index(side_before[1])
                angle = [hexa_before[idx1 - 2], side_before[0], side_before[1], hexa[idx2 + 2]]
                angles.append(angle)
    return angles



impropers = np.vstack((get_angles(hexagons), get_between_hexagons_angles(hexagons, int(columns))))
impropers = impropers + 1






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
for i in range(1, u.atoms.n_atoms+1):
    topology_file.write(f"  {i:<5}     TC5     0     GRA     B{i:<5}     {i:<5}     0     {int(u.atoms[i-1].mass)}\n")


# Bonds

topology_file.write( "\n[ bonds ]\n" )
topology_file.write( "; i	 j	  funct	 length	 kb\n" )
for i in bonds(u):
    topology_file.write(f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.24595     10000\n")










# Angles

topology_file.write( "\n[ angles ]\n" )
topology_file.write( "; i	 j	 k	 funct	 angle	 force_k\n" )
for i in angles(u):
    topology_file.write(f"  {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     1    120     50\n")





# Improper Dihedrals

topology_file.write( "\n[ dihedrals ]\n" )
topology_file.write( "; i	 j	 k	 l     funct	 ref.angle     force_k\n" )
for i in impropers:
    topology_file.write(f"  {i[0]:<3}     {i[1]:<3}     {i[2]:<3}     {i[3]:<3}    180     400\n")
    
    




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




                        
        
        

        


    
        
                

            
        



        


        
