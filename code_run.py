import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
import numpy as np
import argparse


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


u = mda.Universe('input.gro')
c = np.unique(u.atoms.positions[:,1]) # finding the unique y-coordinate so that we can loop through each row of atoms
u.atoms.masses = 36    # Since it's a TC5 bead, we use a mass of 36 for each


'''
Identifying the indices of the virtual sites, and setting their mass to zero
'''
for j in range(len(c)): # Looping through each row
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
    return list_of_angles
                        
        
        

        


    
        
                

            
        



        


        
