
from sec_struct_analysis.sec_struct_analysis_from_vmd import SecondaryStructureAnalysisForVMD as rs


path_to_file = '/Users/aref/work/scratches/test/ssa_test/new/sec_struct_ep2.dat'

file1 = rs(path_to_file)
file1.sec_struct_analysis()
file1.residues_sec_struct()
file1.octave_for_residues(91)
