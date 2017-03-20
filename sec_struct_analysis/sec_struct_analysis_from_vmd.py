import re
import os
import pandas as pd


def sec_struct_analysis(path_to_file):
    """
    analyze the output of the secondary structure calculation in vmd. it
    writes also an input script for octave to plot the quantitative results.
    :param path_to_file: str,
    :return void
    """
    en = path_to_file.split('/')
    file_name = en[-1]
    path = '/'.join(en[:-1])
    os.chdir(path)

    # reading data
    read = []
    with open(file_name, 'r') as inf:
        for l in inf:
            if re.match('^\d+\s', l, re.M):
                read.append(l)
    # preparing the data for analysis
    temp = []
    for l in read:
        line = re.sub('\s', ',', l.strip())
        temp.append(line)
    # empty memory
    del read

    with open('temp.csv', 'w') as ouf:
        for l in temp:
            ouf.write(l + '\n')
    # empty memory
    del temp

    # writing dataframe for the secondary structure per frame
    n = ['resid', 'chain', 'segname', 'frame', 'sec_struct']
    raw_df = pd.read_csv('temp.csv', sep=',', index_col=0, names=n)
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
    for k in dict_sec.keys():
        if len(dict_sec[k]) > 0:
            frame, count = [], []
            # gathering resutls for the secondary structure (k)
            for i in dict_sec[k]:
                frame.append(i[0])
                count.append(i[1].rstrip('\n'))
            # writing these resutls (counts in each frame for the corresponding sec. struct.)
            with open(k + '.csv', 'w') as out:
                for x, y in zip(frame, count):
                    out.write(x + ',' + y + '\n')

            # writing statistical data
            stat.write(k + ":\n")
            df = pd.DataFrame.from_dict(dict_sec[k])
            df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1])
            stat.write("STD of bonds number = {}\n".format(df.iloc[:, 1].std(axis=0).__round__(1)))
            stat.write("Mean of bonds number = {}\n".format(df.iloc[:, 1].mean(axis=0).__round__(1)))
            del df
            # writing octave plotting script
            octave.write(k + "=csvread('" + k + ".csv');\n")
            octave.write("x" + k + "=" + k + "(:,1);\n")
            octave.write("y" + k + "=" + k + "(:,2);\n")

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
