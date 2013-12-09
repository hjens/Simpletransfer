#!/usr/bin/python
import numpy as np
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import pylab as pl
import Simpletransfer1D as st
import sys
from optparse import OptionParser

kpc 	= 3.08567758131e21	#kpc to cm

#Parse command line arguments
parser = OptionParser()
parser.add_option('-p', '--parfile', dest='params_filename', help='Name of parameters file', default='parameters.txt')
parser.add_option('-t', '--trans_file', dest='tdata', help='Name of transmission file', default='transmission.txt')
parser.add_option('-d', '--data_file', dest='data_file', help='Name of data file')
(options,args) = parser.parse_args()

params_filename = options.params_filename
tdata = options.tdata

num_plots = 1 if options.data_file == None else 2
print 'num_plots', num_plots

#Read parameters and transmission data
print 'Reading parameters from:', params_filename
params = st.read_params(params_filename)

#Read redshift
z = float(params['z'])
Hz = 70.*np.sqrt(0.27*(z+1)**3 + (1.-0.27))
print 'Plotting for redshift:', z

data = np.loadtxt(tdata)
wavel = data[:,0]
T = data[:,1]

pl.figure()

#Plot transmission
ax = host_subplot(num_plots,1,1, axes_class = AA.Axes)
ax.plot(wavel, T)
ax.set_xlabel('$\lambda \; \mathrm{[\AA]}$')
ax.set_ylabel('Transmitted fraction')
ax.set_ylim([0,1])

ax2 = ax.twin()

Mpc = np.arange(-50,50,1)
Mpc_lambda = 1215.67*(Mpc*Hz/3.e5+1.)
ax2.set_xticks(Mpc_lambda)
ax2.set_xticklabels(['$'+str(-x)+'$' for x in Mpc])
ax2.set_xlabel('$\mathrm{pMpc}$')

#Plot density data
if num_plots == 2:
	pl.subplot(num_plots, 1, 2)
	col_data = np.loadtxt(options.data_file)

	dist_data = col_data[:,0]/kpc #distance in kpc
	x_hi_data = col_data[:,1]
	temperature_data = col_data[:,2]
	rho_n_data = col_data[:,3] 
	nhi = rho_n_data*x_hi_data

	pl.plot(dist_data, temperature_data/temperature_data.max(), label='$T\; (\\times %.3e \mathrm{K})$' % temperature_data.max())
	pl.plot(dist_data, nhi/nhi.max(), label='$n_{HI}\; (\\times %.3e \mathrm{cm^{-3}})$' % nhi.max())
	pl.legend(loc='best', prop={'size':10})
	pl.xlabel('$d \; \mathrm{kpc}$')

pl.show()
