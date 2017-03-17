import sys
sys.path.append("/home/users/a/aref.takiden/Dropbox/pythonProjects/general_package")
from general_package.new_toolbox import GatheringEnergy
from aMD.acmd_boosting_factors_tools import acmd_tools
import numpy as np


# directory of the production files
directroy = "/work/aref.takiden/acmd/cph1/more-acceleration-from-1run/h470hsp/prod"
ge = GatheringEnergy(directroy)
files = ge.create_lst_of_files(6, 10, "prod", ".", "out")
# grapping the energy lines for the output
ener = ge.gather_energies(files)
# gathering the dihedral and potential energies
dihe = []
potential = []
for l in ener:
    en = l.split()
    dihe.append(float(en[4]))
    potential.append(float(en[13]))
# calculating the average value for the dihedral- and the potential- energies
dihe_aver = np.round(np.average(dihe), 3)
pot_aver = np.round(np.average(potential), 3)
print('======================')
print('Average dihedral energy = ', dihe_aver)
print('Average potential energy = ', pot_aver)
######################################
######################################
# create an object of the class
ac = acmd_tools()

# the directroy to your structure file (*.pdb)
path = "/work/aref.takiden/acmd/cph1/more-acceleration-from-1run/h470hsp/solve_and_neutralize"
struct_file = "namd.pdb"

# number of solute residues
# the third parameter is the segment name, which is by default PROT
resnum = ac.count_residues(path, struct_file, "PROT")
# counting the atoms in the structure
# -> get the names of all segments (1)
# -> get the total number of atoms (2)

# (1)
segnames_list, number_of_segments = ac.get_segment_name_and_number(path, struct_file)

# (2)
counter = 0
for name in segnames_list:
    try:
        counter += ac.count_atoms_of_segment(path, struct_file, name)
    except IOError:
        print('the structure file was not found..')

# geting the factors you want for your ACMD simulation
ac.acmd_boosting_factors(dihe_aver, resnum, pot_aver, counter)


