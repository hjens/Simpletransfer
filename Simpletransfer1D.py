#!/usr/bin/python

#Reads a data file containing distances, densities, temperatures and
#ionized fractions along a line of sight
#Outputs a text file with (wavelength, transmitted fraction)
#Run ./Simpletransfer1D.py -h to see usage 
#Based on IGMTransfer by Peter Laursen

import numpy as np
from scipy.interpolate import interp1d
from optparse import OptionParser


#Physical constants
c		= 2.99792458e10;	#Speed of light; cm/s
nu_0	= 2.46607e15;  		#Lya line-center frequency; Hz
Dnu_L	= 9.936e7			#LyA natural line width; Hz
f_12	= 0.4162          	#Oscillator strength
e    	= 4.803206e-10		#electron charge
m_e		= 9.1093897e-28   	#Electron mass; g
H0  	= 72.#70.				#Hubble constant
Omega0	= 0.044+0.254#0.27				#Matter fraction
lam 	= 1.-Omega0			#Cosmological constant
kpc 	= 3.08567758131e21	#kpc to cm
abu_he	= 0.074				#fractional helium abundance
#--------------

def read_params(params_filename):
	'''Read parameter file and return
	a key/value dictionary '''

	f = open(params_filename, 'r')
	lines = f.readlines()
	f.close()
	params = {}
	for line in lines:
		key, val = line.split('#')[0].split('=')
		params[key.strip()] = val.strip()
	return params

def sigma(Dnu_D, x, idx):
	''' Get LyA cross section at the specified position '''

	a = Dnu_L/(2.*Dnu_D[idx]);

	x2 = x**2
	z = (x2-0.855)/(x2+3.42)
	q = np.where(z <= 0, 0.0, z*(1.+21./x2)*a/(np.pi*(x2+1.0)) * (0.1117 + z*(4.421 + z*(-9.207 + 5.674*z))) )
	voigt = np.sqrt(np.pi)*q + np.exp(-x2)
	sigma_x = f_12 * np.sqrt(np.pi) * e*e/(m_e*c*Dnu_D[idx])*voigt

	return sigma_x

def lorentz_transform(x, Dnu_D, U, idx0):
	oldD = Dnu_D[idx0-1]
	oldU = U[idx0-1]
	newD = Dnu_D[idx0]
	newU = U[idx0]
	
	return (x+oldU)*(oldD/newD) - newU

def init_freq(Dnu_D, params):
	''' Initialize the systemic frequency '''

	wavel = np.linspace(float(params['BWLower']), float(params['BWUpper']), int(params['SpecRes']) )
	nu = c/(wavel * 1.e-8)
	return (nu-nu_0)/Dnu_D[0]

def get_tau(dist, n_HI, Dnu_D, U, params):
	''' Calculate the total tau along the line of sight '''

	#This array will hold the finished tau
	tau = np.zeros(int(params['SpecRes']))

	#Init frequency 
	x = init_freq(Dnu_D, params)

	#Calculate tau
	for i in range(1,len(dist)):
		sigma_x = sigma(Dnu_D, x, i)
		x = lorentz_transform(x, Dnu_D, U, i)
		tau += (dist[i]-dist[i-1])*kpc*n_HI[i]*sigma_x

	return tau

def get_interpolated_array(in_array, new_len, kind='nearest'):
	''' Get a higher-res version of in_array 
	new_len gives the length of the new array, kind gives the
	type of interpolation to perform ''' 

	old_len = len(in_array)
	func = interp1d(np.linspace(0,1,old_len), in_array, kind=kind)
	out_array = func(np.linspace(0,1,new_len))
	return out_array

#-----------------------------
#If you are writing your own script, this is the function you want to use!
def run_sim(params, col_data):
	'''
	Run a LyA transfer simulation
	Parameters:
		* 	params --- dictionary containing simulation parameters. Must include the keys 
			'BWLower', 'BWUpper', 'SpecRes' 'cell_split' and 'z' 
		* col_data --- array with size Nx4. 
			column 0: distance in cm
			column 1: ionized fraction
			column 2: temperature in K
			column 3: total hydrogen density in cm^-3
	Returns:
		tau, wavel
		tau is the optical depth as a function of wavelength
		wavel is the wavelength in Angstrom
	'''
	#Get column data
	dist_data = col_data[:,0]/kpc #distance in kpc
	x_hi_data = col_data[:,1]
	temperature_data = col_data[:,2]
	rho_n_data = col_data[:,3] 

	#Increase the resolution by interpolating the data
	cell_split = int(params['cell_split'])
	new_len = len(dist_data)*cell_split
	dist = get_interpolated_array(dist_data, new_len, kind='linear')
	x_hi = get_interpolated_array(x_hi_data, new_len, kind='nearest')
	temperature = get_interpolated_array(temperature_data, new_len, kind='nearest')
	rho_n = get_interpolated_array(rho_n_data, new_len, kind='nearest')

	Dnu_D = 1.0566e11*(temperature/1.e4)**(1./2.)
	Hz = H0*np.sqrt(Omega0*float(params['z'])**3 + lam)	#Hubble parameter
	U = dist*Hz/1000./12.845/((temperature/1.e4)**(1./2.)) #Gas velocity/thermal width
	n_HI = rho_n*x_hi*(1.-abu_he) #Number density of HI

	tau = get_tau(dist, n_HI, Dnu_D, U, params)
	wavel = np.linspace(float(params['BWLower']), float(params['BWUpper']), int(params['SpecRes']) )

	return tau, wavel

#------ Run as standalone program ---------------


if __name__ == '__main__':
	#Parse command line arguments
	parser = OptionParser()
	parser.add_option('-p', '--parfile', dest='params_filename', help='Name of parameters file', default='parameters.txt')
	parser.add_option('-d', '--data_file', dest='data_file', help='Name of data file')
	(options,args) = parser.parse_args()

	data_filename = options.data_file
	params_filename = options.params_filename

	#Read the data file
	print 'Reading data from ', data_filename
	file_data = np.loadtxt(data_filename)
	print 'Reading parameters from ', params_filename
	params = read_params(params_filename)
	print 'Calculating tau...'

	#Figure out if the data array is in Tom Theun's format
	col_data = np.zeros((file_data.shape[0],4))
	if file_data.shape == col_data.shape:
		col_data = file_data
	else:
		col_data[:,:2] = file_data[:,:2]
		col_data[:,2:] = file_data[:,3:5]


	tau, wavel = run_sim(params, col_data)

	#Save data
	np.savetxt(params['outfile'], np.vstack([wavel, np.exp(-tau)]).transpose())
	print 'Data saved to ', params['outfile']
