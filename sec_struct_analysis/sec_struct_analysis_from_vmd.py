import re
import os
import errno
import collections
import pandas as pd

class SecondaryStructureAnalysisForVMD():
    """
    A toolbox to analyze the output of secondary structure calculation from VMD
    Note:
        -> none of its methods are static, since all related to the same directory, which
        is defined in the __init__ method for each instance!
    """

    def __init__(self, path_to_file):

        en = path_to_file.split('/')
        self.file_name = en[-1]
        self.path = '/'.join(en[:-1])
        try:
            os.chdir(self.path)
        except FileNotFoundError:
            print('\nDirectory:\n{}\nwas not found\n'.format(self.path))
            exit(errno.ENOENT)

        # reading data
        self.read = []
        try:
            with open(self.file_name, 'r') as inf:
                for l in inf:
                    if re.match('^\d+\s', l, re.M):
                        self.read.append(l)
        except FileNotFoundError:
            print('\nThe file\n{}\nwas not found\n'.format(self.file_name))
            exit(errno.ENOENT)
        except IsADirectoryError:
            print('\nPlease give the name of the file, Not its directory\n')
            exit(errno.ENOENT)

    def residues_sec_struct(self):
        """
        analyze the output of the secondary structure calculation in vmd with concentration on
        amino acids. it writes also an input script for octave to plot the quantitative results.
        :param path_to_file:str
        :return: void
        """

        # preparing the data for analysis
        temp = []
        for l in self.read:
            line = re.sub('\s', ',', l.strip())
            temp.append(line)
        # empty memory

        with open('temp.csv', 'w') as ouf:
            for l in temp:
                ouf.write(l + '\n')
        # empty memory
        del temp

        # writing dataframe for the secondary structure per frame
        n = ['resid', 'chain', 'segname', 'frame', 'sec_struct']
        raw_df = pd.read_csv('temp.csv', sep=',', names=n)
        # grouping data
        sec_by_res = raw_df.groupby(['resid'])['sec_struct'].value_counts()
        # empty memory
        os.remove('temp.csv')
        # writing grouped data to a file
        df = sec_by_res.to_frame('counts')
        df['normalized'] = raw_df.groupby(['resid'])['sec_struct'].value_counts(normalize=True).__round__(3)
        del raw_df
        df.to_csv('res_sec_structure.csv')
        return


    def sec_struct_analysis(self):
        """
        analyze the output of the secondary structure calculation in vmd. it
        writes also an input script for octave to plot the quantitative results.
        :param path_to_file: str,
        :return void
        """

        # preparing the data for analysis
        temp = []
        for l in self.read:
            line = re.sub('\s', ',', l.strip())
            temp.append(line)
        # empty memory

        with open('temp.csv', 'w') as ouf:
            for l in temp:
                ouf.write(l + '\n')
        # empty memory
        del temp

        # writing dataframe for the secondary structure per frame
        n = ['resid', 'chain', 'segname', 'frame', 'sec_struct']
        raw_df = pd.read_csv('temp.csv', sep=',', index_col=0, names=n)
        # run statistics on each residues
        # residues_stats(raw_df)
        os.remove('temp.csv')
        # grouping the dataframe by col. frame, then counting the occurence of sec. struct.
        # this will produce a pandas.series.Series
        sec = raw_df.groupby(['frame'])['sec_struct'].value_counts()
        sec.to_csv('sec.csv', sep=',', header=True)
        # getting types of sec. struct.
        sec_set = set(raw_df['sec_struct'])
        # empty memory
        del sec
        del raw_df

        # collecting data
        dict_sec = {a: [] for a in sec_set}
        with open('sec.csv', 'r') as inf:
            for l in inf:
                if re.match('^\d+', l, re.M):
                    en = l.split(',')
                    for k in dict_sec.keys():
                        if en[1] == k:
                            dict_sec[k].append((en[0], en[2]))

        # writing results out (octave script and statistical data)
        # writing script for octave out
        octave = open('oct.m', 'w')
        stat = open('stat.dat', 'w')
        # parameters for octave script
        legend = []
        x_names, y_names = [], []
        sorted_dict_sec = collections.OrderedDict(sorted(dict_sec.items()))
        for k in sorted_dict_sec.keys():
            if len(sorted_dict_sec[k]) > 0:
                frame, count = [], []
                # gathering resutls for the secondary structure (k)
                for i in sorted_dict_sec[k]:
                    frame.append(i[0])
                    count.append(i[1].rstrip('\n'))
                # writing these resutls (counts in each frame for the corresponding sec. struct.)
                with open(k + '.csv', 'w') as out:
                    out.write("Frame, Count\n")
                    for x, y in zip(frame, count):
                        out.write(x + ',' + y + '\n')

                # writing statistical data
                stat.write(k + ":\n")
                df = pd.DataFrame.from_dict(sorted_dict_sec[k])
                df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1])
                stat.write("STD of AA number = {}\n".format(df.iloc[:, 1].std(axis=0).__round__(1)))
                stat.write("Mean of AA number = {}\n".format(df.iloc[:, 1].mean(axis=0).__round__(1)))
                del df
                # writing octave plotting script
                octave.write(k + "=csvread('" + k + ".csv');\n")
                octave.write("x" + k + "=" + k + "(2:end,1);\n")
                octave.write("y" + k + "=" + k + "(2:end,2);\n")

                legend.append(k)
                x_names.append("x" + k)
                y_names.append("y" + k)

        # empty memory for the next data-file
        octave.write('x=[' + ';'.join(x_names) + '];\n')
        octave.write('y=[' + ';'.join(y_names) + '];\n')

        plot = ''
        for i, j in zip(x_names, y_names):
            plot += ','.join([i, j, '"-*",'])

        octave.write('plot(' + plot.rstrip(',') + ')\n')
        octave.write('axis([min(min(x)) max(max(x)) min(min(y))-1 max(max(y))+1])\n')

        line = ''
        for x in legend:
            line += '"' + x + '",'
        octave.write('legend(' + line.rstrip(',') + ')\n')
        octave.write('xlabel("Frame")\n')
        octave.write('ylabel("Number of bonds")\n')
        octave.write('print -dpng -color graph.png')
        # close octave
        octave.close()
        del octave
        del legend
        del x_names
        del y_names

        # close stat-file
        stat.close()

    def octave_for_residues(self, resid):
        """
        extracting the secodary structure data for the given residues from the file res_sec_structure.csv
        :param resid: integer correspond to the consedered residue
        :return: void
        """
        # test the type of argument?
        if not isinstance(resid, int):
            raise ValueError('octave_for_residues method: You have to give an integer..')
            exit(errno.ENOENT)

        # test if input file exist?
        try:
            df = pd.read_csv('res_sec_structure.csv', sep=',', header=0)
        except FileNotFoundError:
            print('\nCould NOT read\nres_sec_structure.csv\nPlease run residues_sec_struct() method'
                  ', then run octave_for_residues() method afterward')
            exit(errno.ENOENT)


        # test if argument valid?
        sdf = df.loc[df['resid'] == resid]
        if len(sdf) != 0:
            print("\nResidue {} data:".format(resid))
            print(sdf)
            print("\n")
            sec_str = list(sdf['sec_struct'])
            perc = list(sdf['normalized'])
            # writing the octave script for plotting
            ouf = open('res_oct.m', 'w')
            line = ''
            for i in perc:
                line += str(i.__round__(3)) + ' '
            ouf.write('y = [' + line + '];\n')
            ouf.write('bar(y)\n')
            xtickes = ""
            for x in sec_str:
                xtickes += "'" + x + "',"
            xtickes = xtickes.rstrip(',')
            ouf.write("set (gca, 'xticklabel', {" + xtickes +"})\n")
            ouf.write('axis ([0,1,0,1]); axis "autox";\n')
            ouf.write('xlabel("Secondary Structure")\n')
            ouf.write('ylabel("Frequency")\n')
            ouf.write('print -dpng -color "residue_sec_struct.png"\n')
            # closing octave script file
            ouf.close()
        else:
            print('\nError: The given resid = {} is not available in res_sec_structure.csv'.format(resid))
            print('Warning: Input code for octave was not written.')
        return

    # def residues_in_sec_structure(self):
