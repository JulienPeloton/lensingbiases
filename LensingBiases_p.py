# Copyright (C) 2016 Lewis, Peloton
########################################################################
# Python script to compute CMB weak lensing biases (N0, N1)
# and derivatives. Internal computations are done in Fortran.
# Authors: (Fortran) Antony Lewis, (Python, and f2py) Julien Peloton
# Contact: j.peloton@sussex.ac.uk
########################################################################
from LensingBiases_f import lensingbiases as lensingbiases_f
from LensingBiases_f import checkproc as checkproc_f

import os
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pylab as pl
pl.ioff()
import argparse

def addargs(parser):
	''' Parse command line arguments '''
	parser.add_argument('-phifile',dest='phifile',
			help='CAMB file containing the fiducial lensing potential',required=True)
	parser.add_argument('-lensedcmbfile',dest='lensedcmbfile',
			help='CAMB file containing the fiducial lensed spectra',required=True)
	parser.add_argument('-FWHM',dest='FWHM',
			help='Beam width (FWHM) in arcmin',required=True,type=float)
	parser.add_argument('-noise_level',dest='noise_level',
			help='Temperature noise level (uk.arcmin). Polar is sqrt(2) bigger.',required=True,type=float)
	parser.add_argument('-lmin',dest='lmin',
			help='Minimum multipole',default=2)
	parser.add_argument('-lmaxout',dest='lmaxout',
			help='Maximum multipole for the output',default=500)
	parser.add_argument('-lmax',dest='lmax',
			help='Maximum multipole for the computation',default=500)
	parser.add_argument('-lmax_TT',dest='lmax_TT',
			help='Maximum multipole for temperature',default=500)
	parser.add_argument('-tmp_output',dest='tmp_output',
			help='Output folder, where files will be written',default='./toto')

def grabargs(args_param=None):
	''' Parse command line arguments '''
	parser = argparse.ArgumentParser(description='Package to compute N0 and N1 lensing biases.')
	addargs(parser)
	args = parser.parse_args(args_param)
	return args

def checkproc_py():
	'''
	Routine to check the number of processors involved
	in the computation (Fortran routines use openmp).
	'''
	nproc = checkproc_f.get_threads()
	if nproc > 1:
		print 'You are using ',nproc,' processors'
	else:
		print '###################################'
		print 'You are using ',nproc,' processor'
		print 'If you want to speed up the computation,'
		print 'set up correctly your number of task.'
		print 'e.g in bash, if you want to use n procs,'
		print 'add this line to your bashrc:'
		print 'export OMP_NUM_THREADS=n'
		print '###################################'

def compute_n0_py(from_args=None,phifile=None,lensedcmbfile=None,
					FWHM=None,noise_level=None,
					lmin=None,lmaxout=None,lmax=None,lmax_TT=None,
					tmp_output=None):
	"""
	Routine to compute the N0 Gaussian bias.
	It calls internally the Fortran routine for speed-up.
	Input:
		* from_args: class, contains all argument for the routine (see addargs).
			If specified, you do not have to specified other arguments.
		* phifile: string, path to the file containing the fiducial lensing potential
		* lensedcmbfile: string, path to the file containing the fiducial lensed CMB spectra
		* FWHM: float, beam width in arcmin
		* noise_level: float, Temperature noise level (uk.arcmin). Polar is sqrt(2) bigger.
		* lmin: int, Minimum multipole
		* lmaxout: int, Maximum multipole for the output
		* lmax: int, Maximum multipole for the computation
		* lmax_TT: int, Maximum multipole for temperature
		* tmp_output: string: Output folder, where files will be written
	Output:
		* return bins, lensing potential, matrix containing
			all N0s, and names of spectra ordered.
	"""
	if from_args is not None:
		lensingbiases_f.compute_n0(args.phifile,args.lensedcmbfile,
			args.FWHM/60.,args.noise_level,
			args.lmin,args.lmaxout,args.lmax,args.lmax_TT,
			args.tmp_output)
	else:
		lensingbiases_f.compute_n0(phifile,lensedcmbfile,
			FWHM/60.,noise_level,
			lmin,lmaxout,lmax,lmax_TT,
			tmp_output)

	n0 = np.loadtxt(os.path.join(args.tmp_output,'N0_analytical.dat')).T
	indices = ['TT','EE','EB','TE','TB','BB']
	bins = n0[0]; phiphi = n0[1]
	n0_mat = np.reshape(n0[2:],(len(indices),len(indices),len(bins)))

	return bins, phiphi, n0_mat, indices

def compute_n1_py(from_args=None,phifile=None,lensedcmbfile=None,
					FWHM=None,noise_level=None,
					lmin=None,lmaxout=None,lmax=None,lmax_TT=None,
					tmp_output=None):
	"""
	Routine to compute the N1 bias.
	It calls internally the Fortran routine for speed-up.
	Input:
		* from_args: class, contains all argument for the routine (see addargs).
			If specified, you do not have to specified other arguments.
		* phifile: string, path to the file containing the fiducial lensing potential
		* lensedcmbfile: string, path to the file containing the fiducial lensed CMB spectra
		* FWHM: float, beam width in arcmin
		* noise_level: float, Temperature noise level (uk.arcmin). Polar is sqrt(2) bigger.
		* lmin: int, Minimum multipole
		* lmaxout: int, Maximum multipole for the output
		* lmax: int, Maximum multipole for the computation
		* lmax_TT: int, Maximum multipole for temperature
		* tmp_output: string: Output folder, where files will be written
	Output:
		* return bins, lensing potential, matrix containing
			all N0s, and names of spectra ordered.
	"""
	if from_args is not None:
		lensingbiases_f.compute_n1(args.phifile,args.lensedcmbfile,
			args.FWHM/60.,args.noise_level,
			args.lmin,args.lmaxout,args.lmax,args.lmax_TT,
			args.tmp_output)
	else:
		lensingbiases_f.compute_n1(phifile,lensedcmbfile,
			FWHM/60.,noise_level,
			lmin,lmaxout,lmax,lmax_TT,
			tmp_output)

	indices = ['TT','EE','EB','TE','TB','BB']
	n1 = np.loadtxt(os.path.join(args.tmp_output,'N1_All_analytical.dat')).T

	bins = n1[0]
	n1_mat = np.reshape(n1[1:],(len(indices),len(indices),len(bins)))

	return bins, n1_mat, indices

def minimum_variance_n0(N0_array,N0_names,checkit=False):
	'''
	Compute the variance of the minimum variance estimator and the associated weights.
	Input:
		* N0_array: ndarray, contain the N0s to combine
		* N0_names: ndarray of string, contain the name of the N0s to combine (['TTTT', 'EEEE', etc.])
	Output:
		* minimum_variance_n0: 1D array, the MV N0
		* weights*minimum_variance_n0: ndarray, the weights for each spectrum (TT, EE, etc.)
		* N0_names_ordered: 1D array, contain the name of the spectra (TT, EE, etc.)
	'''
	N0_array = np.reshape(N0_array,(len(N0_array)**2,len(N0_array[0][0])))
	N0_names_full = ['%s%s'%(i,j) for i in N0_names for j in N0_names]

	## Fill matrix
	sub_vec = [[name,pos] for pos,name in enumerate(N0_names)]
	dic_mat = {'%s%s'%(XY,ZW):[i,j] for XY,i in sub_vec for ZW,j in sub_vec}

	## Build the inverse matrix for each ell
	def build_inv_submatrix(vector_ordered,names_ordered,dic_mat,nsub_element):
		mat = np.zeros((nsub_element,nsub_element))
		for pos,name in enumerate(names_ordered):
			mat[dic_mat[name][0]][dic_mat[name][1]] = mat[dic_mat[name][1]][dic_mat[name][0]] = vector_ordered[pos]
		return np.linalg.pinv(mat)

	inv_submat_array = np.array([build_inv_submatrix(vector_ordered,N0_names_full,dic_mat,len(N0_names)) for vector_ordered in np.transpose(N0_array)])
	inv_N0_array = np.array([ np.sum(submat) for submat in inv_submat_array ])
	minimum_variance_n0 = 1./inv_N0_array

	weights = np.array([ [np.sum(submat[i]) for submat in inv_submat_array] for i in range(6) ])

	if checkit:
		print 'Sum of weights = ', np.sum(weights*minimum_variance_n0)/len(minimum_variance_n0)
		print 'Is sum of weights 1? ',np.sum(weights*minimum_variance_n0)/len(minimum_variance_n0) == 1.0

	return minimum_variance_n0, weights*minimum_variance_n0

def minimum_variance_n1(bins,N1_array,weights_for_MV,spectra_names,bin_function=None):
	'''
	Takes all N1 and form the mimimum variance estimator.
	Assumes N1 structure is coming from Biases_n1mat.f90
	Input:
		* N1: ndarray, contain the N1 (output of Biases_n1mat.f90)
		* weights_for_MV: ndarray, contain the weights used for MV
		* spectra_names: ndarray of string, contain the name of the spectra ordered
	'''
	## Ordering: i_TT=0,i_EE=1,i_EB=2,i_TE=3,i_TB=4, i_BB=5 (from Frotran)
	names_N1 = ['%s%s'%(i,j) for i in spectra_names for j in spectra_names]

	if bin_function is not None:
		n1_tot = np.zeros_like(bin_centers)
	else:
		n1_tot = np.zeros_like(weights_for_MV[0])

	for estimator_name in names_N1:
		## Indices for arrays
		index_x = spectra_names.index(estimator_name[0:2])
		index_y = spectra_names.index(estimator_name[2:])

		## Interpolate N1 if necessary
		n1_not_interp = N1_array[index_x][index_y]
		if bin_function is not None:
			n1_interp = np.interp(bin_centers,bins,n1_not_interp)
		else:
			n1_interp = n1_not_interp

		## Weights
		wXY_index = spectra_names.index(estimator_name[0:2])
		wZW_index = spectra_names.index(estimator_name[2:4])

		## Update N1
		if bin_function is not None:
			n1_tot += bin_function(weights_for_MV[wXY_index]) * bin_function(weights_for_MV[wZW_index]) * n1_interp
		else:
			n1_tot += weights_for_MV[wXY_index] * weights_for_MV[wZW_index] * n1_interp
	return n1_tot

def plot_biases(bins,phiphi,MV_n0,MV_n1=None,N0_array=None,N1_array=None):
	'''
	Quick plot for inspection
	'''
	tphi = lambda l: (l+0.5)**4/(2.*np.pi) # scaling to apply to cl_phiphi when plotting.
	colors = lambda i: matplotlib.cm.jet(i*60)

	## Plot lensing
	pl.loglog(bins,phiphi,color='grey',label='Lensing')

	## Plot N0
	pl.loglog(bins,MV_n0*tphi(bins),color='black',lw=2,label='N0 bias')
	if N0_array is not None:
		indices = ['TT','EE','EB','TE','TB','BB']
		for i in range(len(N0_array)):
			pl.loglog(bins,N0_array[i][i][:]*tphi(bins),color=colors(i),
					lw=2,alpha=0.2,label=indices[i]+indices[i])

	## Plot N1
	if MV_n1 is not None:
		pl.loglog(bins,MV_n1*tphi(bins),color='black',lw=2,ls='--',label='N1 bias')
	if N1_array is not None:
		indices = ['TT','EE','EB','TE','TB','BB']
		for i in range(len(N1_array)):
			pl.loglog(bins,N1_array[i][i][:]*tphi(bins),color=colors(i),
					ls='--',lw=2,alpha=0.2,label=indices[i]+indices[i])

	pl.xlabel('$\ell$',fontsize=20)
	pl.ylabel(r"$[\ell(\ell+1)]^2/(2\pi)C_\ell^{\phi^{XY} \phi^{ZW}}$",fontsize=20)
	leg=pl.legend(loc='best',ncol=2,fontsize=12.5)
	leg.get_frame().set_alpha(0.0)
	pl.savefig('Biases.pdf')
	pl.clf()

if __name__ == "__main__":
	args_param = None
	args = grabargs(args_param)

	## Check openmp
	checkproc_py()

	## Compute N0s, and form MV
	## Example with argparse
	bins, phiphi, n0_mat, indices = compute_n0_py(from_args=args)

	## Example with direct arguments
	# bins, phiphi, n0_mat, indices = compute_n0_py(from_args=None,phifile=args.phifile,lensedcmbfile=args.lensedcmbfile,
	# 					FWHM=args.FWHM,noise_level=args.noise_level,
	# 					lmin=args.lmin,lmaxout=args.lmaxout,lmax=args.lmax,lmax_TT=args.lmax_TT,
	# 					tmp_output=args.tmp_output)
	MV_n0, weights = minimum_variance_n0(n0_mat,indices,checkit=False)

	## Compute N0s, and form MV
	bins, n1_mat, indices = compute_n1_py(from_args=args)
	MV_n1 = minimum_variance_n1(bins,n1_mat,weights,indices,bin_function=None)

	plot_biases(bins,phiphi,MV_n0,MV_n1=MV_n1,N0_array=n0_mat,N1_array=n1_mat)
