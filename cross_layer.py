import numpy as np
import os
import MDAnalysis as mda
import pandas as pd

def get_zpositions(u):
    """put the selected atoms z positions in a pandas dataframe"""
    # define a dict to store the z positons and their atomic index
    atom_positions = {}
    # iterate over the traj
    for ts in u.trajectory:
        # iterate over all F ions
        for atomindex in range(36, 108):
            tempz = ts.positions[atomindex][2]
            # if index in the key, append the position
            # else, create another key:value pair and put the vlaue as list
            if atomindex in atom_positions:
                atom_positions[atomindex].append(tempz)
            else:
                atom_positions[atomindex] = [tempz]

    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(atom_positions)
    return df

def count_one_atom_cross(atomzposlist):
    """this function counts how many times one atom cross Sn-Sn layer
    across the trajectory
    atomzposlist: z positions of one atom
    atomindex: F ion index to calculate
    """
    # initiate a variable to accumulate the crossing number
    atomcross = 0
    # count the crossing for one atom
    # we define 0 and 1 as threshold for crossings
    for i in range(1, len(atomzposlist)):
        if atomzposlist[i] > 0 and atomzposlist[i - 1] < 0:
            atomcross += 1
        elif atomzposlist[i] < 0 and atomzposlist[i - 1] > 0:
            atomcross += 1
        elif atomzposlist[i] < 10.9 and atomzposlist[i - 1] > 10.9:
            atomcross += 1
        elif atomzposlist[i] > 10.9 and atomzposlist[i - 1] < 10.9:
            atomcross += 1
    return atomcross

path2dir = os.getcwd()
udata = mda.Universe(os.path.join(path2dir, 'XDATCAR.pdb'))
zposdf = get_zpositions(udata)
crossings = []
for column_name in zposdf.columns:
    atomcross = count_one_atom_cross(zposdf[column_name])
    crossings.append(atomcross)

print('average crossing times are', sum(crossings)/len(crossings))


# use the index when atom cross the layer to define the duration
def get_one_atom_continuous_stay(atomzposlist):
    """count for how long the ion will stay at one layer before moving to the
    next layer
    """
    durations = []
    duration_start_index = 0
    for i in range(1, len(atomzposlist)):
        if atomzposlist[i] > 0 and atomzposlist[i - 1] < 0:
            durations.append(i - duration_start_index)
            duration_start_index = i
        elif atomzposlist[i] < 0 and atomzposlist[i - 1] > 0:
            durations.append(i - duration_start_index)
            duration_start_index = i
        elif atomzposlist[i] < 10.9 and atomzposlist[i - 1] > 10.9:
            durations.append(i - duration_start_index)
            duration_start_index = i
        elif atomzposlist[i] > 10.9 and atomzposlist[i - 1] < 10.9:
            durations.append(i - duration_start_index)
            duration_start_index = i
    return durations

# calculate the stay for all atoms
all_atoms_stay = []
for column_name in zposdf.columns:
    all_atoms_stay.extend(get_one_atom_continuous_stay(zposdf[column_name]))
print('the average staying steps are', sum(all_atoms_stay)/len(all_atoms_stay))
