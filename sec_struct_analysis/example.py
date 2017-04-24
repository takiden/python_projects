import sys
sys.path.append("/Users/aref/python_projects")
from sec_struct_analysis.sec_struct_analysis_from_vmd import SecondaryStructureAnalysisForVMD as rs

# directory to the data file from VMD
path_to_file='/Users/aref/work/agp1/exostructM15/mutants/wc_correct/analysisData_w2c/sec-struct/w2c-secstruct-438to474.dat'
# starting the class as file1
file1 = rs(path_to_file)
# run a secondary structure analysis
file1.sec_struct_analysis()
# run an analysis for each residue
file1.residues_sec_struct()
# generate octave code to plot the results
#file1.octave_for_residues(91)
