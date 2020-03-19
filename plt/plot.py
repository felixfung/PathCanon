#!/usr/bin/python

##############################################################################
# plot.py
#
# plot one or more, or all timeseries
# given a text file with format:
# each column, white space separated, is a timeseries
# first line is header for the timeseries names
# and the first column is "Time"
# all other entries are numerical
#
# usage:
#
# ./plot.py timeseries.txt fieldname1 fieldname2 ...
#
# fieldname can be "All" to print all fields
#
##############################################################################

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]

with open( filename ) as file:
    header = file.readline()
fields = header.split()

fieldname2idx = {}
for i in range(len(fields)):
    fieldname2idx[ fields[i] ] = i 

data = np.loadtxt( filename, skiprows=1 )

print_all = False
plot_field = { fieldname: False for fieldname in fields }
for i in range(2,len(sys.argv)):
    if sys.argv[i] == 'All':
        print_all = True
    plot_field[ sys.argv[i] ] = True

for field in fields:
    if print_all or plot_field[field]:
        if field == 'Time':
            continue
        plt.figure()
        plt.plot( data[:,0], data[:,fieldname2idx[field]] )
        plt.xlabel('Time')
        plt.ylabel(field)
        plt.xlim([ data[:,0].min(), data[:,0].max() ])
        plt.ylim([ data[:,fieldname2idx[field]].min(), \
                   data[:,fieldname2idx[field]].max() ])
plt.show()
