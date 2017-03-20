import re, os, sys


class GatheringEnergy:
    """
    this very short class is simply to create a list of files that are in a directory
     and then extract the energy terms out of these files
    """
    @staticmethod
    def gather_energies(path, lst_of_files):
        """
        gathered the energies from the list of files and write them into one file
        and return void
        :return list
        """
        try:
            os.chdir(path)
        except FileNotFoundError:
            print('\nDirectory \n\t{}\nwas not found.\nTerminating\n'.format(path))
            exit(1)

        energies = []
        damaged_files = {}
        for f in lst_of_files:
            try:
                damaged_files[f] = []
                with open(f, "r") as inf:
                    for l in inf:
                        mat = re.match(r'^ENERGY:', l, re.M)
                        if mat:
                            # check if the matched line has the correct entries.
                            en = l.split()
                            for x in en[1:]:
                                try:
                                    float(x)
                                except ValueError:
                                    damaged_files[f].append(l)
                                    break
                            else:
                                energies.append(l)
            except IOError:
                print('file {} was not found..\nTerminating.'.format(f))
                exit(1)
        # check if there are some damaged lines
        damaged = False
        damaged_names = []
        for k, v in damaged_files.items():
            if len(v) > 0:
                damaged = True
                print('\nenergy lines could not be parsed correctly in file {}:\n'.format(k))
                n = k.split('.')[0]+'.corrupted'
                damaged_names.append(n)
                # writing out the damaged lines to the corresponding files
                with open(n, 'w') as ouf:
                    for x in v:
                        ouf.write(x+'\n')

        if damaged:
            print('\nTo see the corrupted files, please type ls *.corrupted \n')
            exit(1)
        else:
            return energies

    @staticmethod
    def renumber_timesteps_of_energies(path, lst_of_energies, config_name):
        """
        reads the energies form function: gather_energies and the configuration file
        :param path: str
        :param lst_of_energies: list
        :param config_name: string
        :return: list
        """
        try:
            os.chdir(path)
        except FileNotFoundError:
            print('\nDirectory \n\t{}\nwas not found.\nTerminating\n'.format(path))
            exit(1)

        try:
            inf = open(config_name, "r")
            read = inf.readlines()
            inf.close()
        except IOError:
            print("confugration file was not found...")

        # setting up the needed variables to renumber the energies
        freq = 0.0
        timestep = 0.0
        for l in read:
            if re.match(r'^outputEnergies', l, re.I):
                en = l.split()
                freq = float(en[1])
            elif re.match(r'^timestep', l, re.I):
                en = l.split()
                timestep = float(en[1])
            elif timestep != 0 and freq != 0:
                break

        t = 0.0

        final = []
        for el in lst_of_energies:
            line = re.sub(r'^ENERGY:\s+\d+\s+', '', el)
            t += freq * timestep * 1e-6
            final.append(str(t.__round__(2)).ljust(6) + line)

        return final

    @staticmethod
    def create_lst_of_files(path, start_num_of_files, last_num_of_files, prefix, seperator, extension):
        """
        create a list of arbitrary number of files that exists in a folder
        :param path: str
        :param start_num_of_files: int
        :param last_num_of_files: int
        :param prefix: string
        :param seperator: string
        :param extension: string
        :return: list
        """
        try:
            os.chdir(path)
            print(path)
        except FileNotFoundError:
            print('\nDirectory \n\t{}\nwas not found.\nTerminating\n'.format(path))
            exit()

        files = []
        i = start_num_of_files
        if isinstance(last_num_of_files, int):
            # for i in range(1,numOfFiles,1):
            while i <= last_num_of_files:
                f = prefix + str(i) + seperator + extension
                if os.path.exists(f):
                    files.append(f)
                else:
                    print('\nFile {} does not exist\nTerminating\n'.format(f))
                    exit(1)
                i += 1
        elif len(last_num_of_files) == 0:
            files.append(prefix + seperator + extension)
        return files


class PdbTools:
    """
    this class contains tools that are based on a pdb files, and for this case, that
    directory is a common feature for all functions, for this reason I will initiate the
    class with the directory itself
    """

    def __init__(self, path, file_name=None):
        if file_name is not None:
            self.path = path
            self.file_name = file_name
            try:
                with open(os.path.join(self.path, self.file_name)) as inf:
                    read = inf.readlines()
                    inf.close()
                    self.lst = read
            except IOError:
                print('the {} was not found'.format(self.file_name))
                exit()
        else:
            self.path = path
            self.file_name = None
            try:
                os.path.exists(self.path)
            except IOError:
                print('the path is not valid')

    def get_path(self):
        return str(self.path)

    def get_file_name(self):
        return str(self.file_name)

    def count_residues(self, file_name=None):
        """
        count the residue in a pdb written by NAMD
        the default segment name ist "PROT" -- Case sensitive --
        :param file_name: string
        :return: integer; the number of residues
        """
        os.chdir(self.path)
        read = []
        try:
            try:
                self.file_name is not None
                inf = open(self.file_name, "r")
                for l in inf:
                    if re.search(r'\sPROT', l, re.M):
                        read.append(l)
                inf.close()
            except IOError:
                file_name = input('please enter the name of the file:\n')
                inf = open(file_name, "r")
                for l in inf:
                    if re.search(r'\sPROT', l, re.M):
                        read.append(l)
                inf.close()
        except IOError:
            print('\nFile was not found\n')
            exit(1)

        if len(read) == 0:
            print('\nCould not parse the segment name -_-\n')
            exit(1)

        # gather lines of the same residues
        residue = []
        for l in read:
            residue.append(l[17:27])

        set_residues = set(residue)
        return len(set_residues)

    def count_atoms_of_segment(self, file_name, segname):
        """
        Count the number of atoms in a pdb file written by NAMD
        The counting based on the segment name
        :param segname:
        :param file_name:
        :return: integer; the number of atoms
        """
        os.chdir(self.path)
        read = []
        try:
            inf = open(file_name, "r")
            pattern = "\s" + segname
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

    def get_segment_name_and_number(self, file_name):
        """
        to get the number of segments in a pdb file and their names
        :param file_name:
        :return: typle(list of segments name, total number of segments )
        """
        os.chdir(self.path)
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

    def get_seq_from_pdb(self):
        """
        to get the sequence of residues in a pdb files, where the record lines
        start with either ATOM or HETATM
        :param file_name: self.file_name
        :return: list of the sequence
        """
        # output list
        seq = []
        # reading list
        read = []
        # reading and processing
        try:
            inf = open(os.path.join(self.path, self.file_name), "r")
            for l in inf:
                if re.match(r'^ATOM\s+', l, re.M):
                    read.append(l[17:26])
                elif re.match(r'^HETATM\s+', l, re.M):
                    read.append(l[17:26])
        except IOError:
            print("\nFile was not found\n")
            exit(1)

        el = read[0]
        entry = read[0].split()[0]
        seq.append(entry)
        i = 1
        while i < len(read):
            if el == read[i]:
                # just go a step further
                i += 1
            elif el != read[i]:
                # get the new entry
                entry = read[i].split()[0]
                seq.append(entry)
                # change the old value of el
                el = read[i]
                i += 1

        return seq

    def one_letter_seq_format(self, lst):
        """
        takes a list of 3-letters sequence and convert it to one-letter sequnce
        :param lst: list.
        :return:  string, string of one-letter sequence
        """

        trans_dic = {"ARG": "R", "HIS": "H", "LYS": "K", "ASP": "D", "GLU": "E",
                     "SER": "S", "THR": "T", "ASN": "N", "GLN": "Q", "CYS": "C",
                     "SEC": "U", "GLY": "G", "PRO": "P", "ALA": "A", "VAL": "V",
                     "ILE": "I", "LEU": "L", "MET": "M", "PHE": "F", "TYR": "Y",
                     "TRP": "W", "HSD": "H", "HSE": "H", "HSP": "H"}

        formated_seq = ""
        for en in lst:
            formated_seq += trans_dic[en]

        return formated_seq

    def get_resid(self, line=None):
        """
        just get the resid out of a line or from the "self.lst variable"
        :param line:str (optional)
        :return:str/list
        """
        if line is not None:
            resid = line[22:26].strip()
            return resid
        else:
            try:
                residlst = []
                # check if a file is given by initialization or not
                len(self.lst) > 0
                for l in self.lst:
                    residlst.append(l[22:26])
                return residlst
            except IOError:
                print('wrong argument. Either you name your file by initializing the instance, or you '
                      'provide a line from a pdb-file')

    def get_segment(self, segname):
        """
        :param segname: string: the name of the segment you want
        :param lst: a list to read the data from
        :return: list: constaining the segment lines
        """
        segment = []
        for l in self.lst:
            if l[72:76].strip() == segname:
                segment.append(l)
        return segment

    def get_residue(self, resid, segname):
        """
        get all corresponding lines of a residue from the pdb file
        :param resid: int
        :param segname: lst
        :return: lst
        """
        id_ = str(resid)
        residue = []
        for l in segname:
            if id_ == self.get_resid(l):
                residue.append(l)

        return residue

    def get_index(self, atom_name):
        """

        :param lst: a list to read the data from
        :param atom_name: searched atom's name
        :return: int: the pdb index of the atom name
        """
        for l in self.lst:
            if l[12:16].strip() == atom_name:
                return l[6:11].strip()
        return 'Index was not found'

    def atom_info(self, atom_name, resid, segname):
        """
        get the line of the corresponding atom in the pdb file
        :param atom_name: str
        :param resid: int
        :param segname: str
        :return: str
        """
        seg_block = self.get_segment(segname)
        residue = self.get_residue(resid, seg_block)

        for l in residue:
            if l[12:16].strip() == atom_name:
                return l
        return ""

    def atom_dict(self, line):
        """
        A line from pdb fle that is generated from NAMD or CHARMM is expected.
        Line format will NOT be checked!
        :param line: str, line from pdb file
        :return: dict,
            dict keys are: index int, atom_name str, resname str, resid int, x float, y float,
            z float, segname str
        """
        dict_infos = {
            'index': int(line[6:11].strip()),
            'atom_name': line[12:16].strip(),
            'resname': line[17:21].strip(),
            'resid': int(line[22:26].strip()),
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54]),
            'segname': line[72:76].strip()
        }
        return dict_infos

    def get_atom_line_by_name(self, line, atom_name):
        """
        read the lines of the pdb-file and return the matching name
        :param line: str,
        :param atom_name: str,
        :return: str, the line containing the atom_name
        """
        if atom_name == line[12:16].strip():
            return line

    @staticmethod
    def get_coor(atom_info_line):
        """
        :param atom_info_line: the pdb-line of the atom
        :return: list of floats, containing the x, y and z coordinates of the selected atoms
        """
        line = atom_info_line
        coor = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        return coor

    @staticmethod
    def calculate_angle(coor_lst1, coor_lst2, coor_lst3):
        """
        you need the coordinates of the three atoms
        :param coor_lst1: lst, list of floats
        :param coor_lst2: lst, list of floats
        :param coor_lst3: lst, list of floats
        :return: float, angle in degrees
        """
        import numpy as np

        a = np.array(coor_lst1)
        b = np.array(coor_lst2)
        c = np.array(coor_lst3)

        ba = a - b
        bc = c - b

        consine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(consine_angle)

        return np.degrees(angle).round(2)

    def calculate_distance(self, atom1, atom2):
        """
        calculate distance between two point in cartesian coordinates in 3D space
        :param atom1: string: line of a pdb file
        :param atom2: string: line of a pdb file
        :return: float: distance in angestrom
        """
        # TODO treate the empty string form atom_info method

        # coordinates of the first atom
        x1 = float(atom1[30:38].strip())
        y1 = float(atom1[38:46].strip())
        z1 = float(atom1[46:54].strip())
        # coordinates of the second atom
        x2 = float(atom2[30:38].strip())
        y2 = float(atom2[38:46].strip())
        z2 = float(atom2[46:54].strip())

        import math
        dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
        return dist.__round__(3)

    def calculate_distance_between_rings_pyrrole_water(self, file_name):
        """
        this function calculate the distance between the pyrrole water and the nitrogens of the rings A, B and C
        of the biliverdin cofactor
        :rtype: void
        :param file_name: string
        :return: void, it prints the values out
        """
        os.chdir(self.path)
        read_opt = []
        inf = open(file_name, "r")
        read_opt = inf.readlines()
        inf.close()
        cw_opt = self.get_segment("CW", read_opt)
        bv_opt = self.get_segment("BV", read_opt)

        # ring B
        n_d = self.atom_info("N_D", bv_opt)
        # ring A
        n_c = self.atom_info("N_C", bv_opt)
        # ring C
        n_a = self.atom_info("N_A", bv_opt)

        cw112 = self.get_residue_by_id(112, cw_opt)
        oh2 = self.atom_info("OH2", cw112)

        print('N(ringA) --> OH2 = ', self.calculate_distance(n_c, oh2))
        print('N(ringB) --> OH2 = ', self.calculate_distance(n_d, oh2))
        print('N(ringC) --> OH2 = ', self.calculate_distance(n_a, oh2))
        return

    def cal_dihedral(self, p):
        """Praxeolitic formula
        1 sqrt, 1 cross product
        this funciton is copied from stackoverflow:http://stackoverflow.com/a/34245697/2804070

        p should be a list of np.array ==> p = [np.array, np.array, ..., etc]
        """
        import numpy as np

        p0 = p[0]
        p1 = p[1]
        p2 = p[2]
        p3 = p[3]

        b0 = -1.0 * (p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2

        # normalize b1 so that it does not influence magnitude of vector
        # rejections that come next
        b1 /= np.linalg.norm(b1)

        # vector rejections
        # v = projection of b0 onto plane perpendicular to b1
        #   = b0 minus component that aligns with b1
        # w = projection of b2 onto plane perpendicular to b1
        #   = b2 minus component that aligns with b1
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1

        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        return np.degrees(np.arctan2(y, x))

class StaticTools():
    """
    collection of usable functions
    """
    import os
    import re
    @staticmethod
    def dssp_output_analyser(path, file_name):

        os.chdir(path)

        seq = []
        dssp = []
        # reading an input file
        try:
            with open(file_name, 'r') as inf:
                for l in inf:
                    en = l.split()
                    seq.append(en[0])
                    if len(en)> 1:
                        jj = en[-1].strip()
                        if jj.isalpha():
                            dssp.append(jj)
        except IOError:
            print('the file ( {} ) was not found..'.format(file_name))
        # processing the input data
        aa = ""
        for l in seq:
            for c in l:
                if c.isalpha():
                    aa += c
        print(aa)
        print(len(aa))
        sec = ""
        for l in dssp:
            for c in l:
                if c.isalpha():
                    sec += c

        # counting the occurrence of each secondary structure
        occu_dict = {}
        for s in set(sec):
            occu_dict[s] = [sec.count(s)]
            print('count of {} = {}'.format(s, sec.count(s)))
        # counting percentage of each sec. structure
        total = 0
        for k, v in occu_dict.items():
            perc = ((v[0] / len(aa)) * 100).__round__(0)
            occu_dict[k].append(perc)
            total += perc
            print('{} = {}%'.format(k, occu_dict[k][1]))
        print('Percentage of defined secondary structure = {}%'.format(total))
        return

    @staticmethod
    def swap(self, a, b):
        return b, a

    @staticmethod
    def convert_whitespace2csv(path, file_name):
        """
        convert the whitespace speperated columns to csv
        :param path: str
        :param file_name:str
        :return: list (csv records)
        """

        os.chdir(path)
        out = []
        with open(file_name, 'r') as inf:
            for l in inf:
                en = l.split()
                out.append(','.join(en))

        return out

