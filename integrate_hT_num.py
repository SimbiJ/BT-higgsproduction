#!/usr/bin/env python3
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import numpy as np
import lhapdf
from scipy import integrate


R_DTYPE = np.float64
C_DTYPE = np.complex128

START      = 1000
NUM_POINTS = 300
STEP       = 1000

# in GeV
m_H    = 125
invmsq = 0.0000001
v      = 250
S      = (13.6*10**3)**2

# load pdf from lhapdf and integral library
GLUON_ID = 21
PDF = lhapdf.mkPDF("PDF4LHC15_nnlo_30_pdfas", 0)

HEAVY_TOP = IntegralLibrary('heavyTop/heavyTop_pylink.so')



def lumi(x1, s12, mu):
	return PDF.xfxQ2(GLUON_ID, x1, mu**2)/x1*PDF.xfxQ2(GLUON_ID, s12/(x1*S), mu**2)/s12



def get_integral_coeff_O0(real_params, complex_params=None):
	if complex_params is None:
		complex_params = []

	# readout results of integral for given parameters
	str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = HEAVY_TOP(
		real_parameters    = real_params,
		complex_parameters = complex_params,
		verbose=True
	)

	# convert results
	integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
	integral_with_prefactor = sp.sympify(integral_with_prefactor.replace('+/-','*value+error*'))
	coeff_O0 = integral_with_prefactor.coeff('eps',0).coeff('value')

	return np.real_if_close(coeff_O0)



def sigma(s12, tau, alpha, C_0):
	coeff = np.pi*alpha**2/(16*(4*np.pi)**2*v**2*s12)
	F_12  = np.abs(2/invmsq*(2+(4/invmsq-s12)*C_0))**2

	integrated_lumi = np.array(
		[integrate.quad(lumi, t, 1., args=(s, m_H))[0] for s, t in zip(s12, tau)],
		dtype=R_DTYPE
	)

	return coeff*F_12*integrated_lumi



def generate_data():
	num_real_param = int(HEAVY_TOP.info['number_of_real_parameters'])
	num_complex_param = int(HEAVY_TOP.info['number_of_complex_parameters'])

	s12 = START + STEP*np.arange(NUM_POINTS, dtype=R_DTYPE)

	# generate list of real parameters
	real_points = np.empty((NUM_POINTS, num_real_param), dtype=R_DTYPE)
	# write data in columns
	real_points[ :, 0] = s12
	real_points[ :, 1] = invmsq

	result = np.array(list(map(get_integral_coeff_O0, real_points)), dtype=C_DTYPE).transpose()

	return s12, np.real_if_close(result)


	
def main():
	s12, result = generate_data()

	result = np.array(result, dtype=R_DTYPE)
	alpha  = PDF.alphasQ2(m_H**2)
	tau    = s12/S

	sigma_hT = sigma(s12, tau, alpha, result)

	s12.tofile('sigma_num_hT_s12.bytes')
	sigma_hT.tofile('sigma_num_hT.bytes')



if __name__ == "__main__":
	main()
