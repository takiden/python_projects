import os
import re


class acmd_tools:
    """
    I use this class to calculate the boosting factors for the accelerated molecular
    dynamic
    """
    def acmd_boosting_factors(self, dihe_aver_ener, solute_res_num, pot_aver_ener, atoms_num):
        """
        calculating the boosting factors for an ACMD simulation
        :param dihe_aver_ener: float
        :param resnum: int
        :param pot_aver_ener: float
        :param atoms_num: int
        :return: void, it prints the results in the terminal
        """

        # explicit the type of the variables
        resnum = int(solute_res_num)
        dihe_aver = float(dihe_aver_ener)
        # dihedral
        # energy boost
        dihe_energy_bosst = (dihe_aver + (4 * resnum)).__round__(3)
        # alpha boost
        dihe_alpha = (1 / 5 * (4 * resnum)).__round__(3)

        # total potential
        # explicit the type of the variables
        a_num = int(atoms_num)
        pot_aver = float(pot_aver_ener)
        # energy boost
        # pot_energy_bosst = (pot_aver + (0.16 * a_num)).__round__(3)
        pot_energy_bosst = (pot_aver + (0.3 * a_num)).__round__(3)
        # alpha boost
        # pot_alpha = (0.16 * atoms_num)
        pot_alpha = (0.1 * atoms_num)
        # print the results in the terminal
        print('===============================')
        print('Dihedral boost\n')
        print('Energy = ', dihe_energy_bosst)
        print('alpha = ', dihe_alpha)
        print('===============================')
        print('Potential energy boost\n')
        print('Energy = ', pot_energy_bosst)
        print('alpha = ', pot_alpha)
        print('===============================')

    def count_residues(self, path, file_name, segment):
        """
        count the residue in a pdb written by NAMD
        the default segment name ist "PROT" -- Case sensitive --
        :param path: string
        :param file_name: string
        :return: integer; the number of residues
        """
        os.chdir(path)
        if segment is None:
            segname = "PROT"
        else:
            segname = segment

        read = []
        pattern = "\s" + segname
        try:
            inf = open(file_name, "r")
            for l in inf:
                if re.search(pattern, l, re.M):
                    read.append(l)
            inf.close()
        except IOError:
            print('\nFile' + file_name + ' was not found\n')
            exit(1)

        if len(read) == 0:
            print('\nCould not parse the PROTEIN segment name -_-\n')
            exit(1)

        # gather lines of the same residues
        residue = []
        for l in read:
            residue.append(l[17:27])

        set_residues = set(residue)
        return len(set_residues)

    def count_atoms_of_segment(self, path, file_name, seg):
        """
        Count the number of atoms in a pdb file written by NAMD
        The counting based on the segment name
        :param path:
        :param file_name:
        :return: integer; the number of atoms in the givin segment
        """
        os.chdir(path)
            
        read = []
        try:
            inf = open(file_name, "r")
            pattern = "\s" + seg
            for l in inf:
                if re.search(pattern, l, re.M):
                    read.append(l)
            inf.close()
        except IOError:
            print('\nFile was not found\n')
            exit(1)

        if len(read) == 0:
            print('\nCould not parse the segment name -_-\n')
            exit(1)
        else:
            return len(read)

    def get_segment_name_and_number(self, path, file_name):
        """
        to get the number of segments in a pdb file and their names
        :param path:
        :param file_name:
        :return: typle(list of segments' names, total number of segments )
        """
        os.chdir(path)
        read = []
        try:
            inf = open(file_name, "r")
            for l in inf:
                if re.match(r'^ATOM\s', l, re.M):
                    read.append(l)
                elif re.match(r'^HETATM\s', l, re.M):
                    read.append(l)
            inf.close()
        except IOError:
            print('\nFile was not found\n')
            exit(1)

        segname = []
        for l in read:
            segname.append(l[72:76])

        set_segment_name = set(segname)
        number_of_segments = len(set_segment_name)
        return list(set_segment_name), number_of_segments
