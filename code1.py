import MDAnalysis as mda
import numpy as np
import nglview as nv
import argparse
positions = []
dist = 0

while True:
    rows = int(input("Enter the number of rows, make sure that it is odd, and greater than 3: "))
    if (rows >=3) and (rows % 2 !=0):
    # the correct input was entered, so we break out of the loop
        break
    else:
    # the incorrect input was entered, so we continue the loop
        continue

while True:
    columns = int(input("Enter the number of beads across a row, make sure that it is odd: "))
    if columns %2 !=0:
        break
    else:
        continue
        



for i in range(rows):
    if i % 2 == 0:
        
        if (i == rows - 1) or (i == 0):
            x = np.arange(0, columns*2.56, 2.56)
            d_array = np.array(np.arange(3, columns + 1, 3) - 1)
            x = np.delete(x,  d_array)
        else:
           x = np.arange(0, columns * 2.56, 2.56)
           
        for k in range(len(x)):
            positions.append([x[k], dist, 0])
        dist+=2.217
    else:
        x = np.arange(-1.28, columns*2.56, 2.56)
        for k in range(len(x)):
            positions.append([x[k], dist, 0])       
        
        dist+=2.217
        

        

    
w = mda.Universe.empty(n_atoms = len(positions), trajectory = True)
w.atoms.positions = positions
w.atoms.write('input.gro')



        
