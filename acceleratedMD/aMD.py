import sys

sys.path.append("..")  # adding python_projects folder to PYTHONLIBRARY
import numpy as np
from acceleratedMD.acmd_boosting_factors_tools import acmd_tools
from general.new_toolbox import GatheringEnergy as ge

# directory to the production files
directory = "/Users/aref/work/scratches/test/prod"

# directory to the starting structure (e.g.: namd.pdb)
path = "/Users/aref/work/scratches/test/struct"

prot_segname = ""

# select the steps you want to calculate their corresponding average energy
files = ge.create_lst_of_files(directory, 6, 10, "prod", ".", "out")

########################################
######## Code               ############
########################################
# grapping the energy lines for the output
ener = ge.gather_energies(directory, files)
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
struct_file = "namd.pdb"

# number of solute residues
# the third parameter is the segment name, which is by default PROT
resnum = ac.count_residues(path, struct_file, prot_segname)
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
