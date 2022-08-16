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
m_H = 125
v   = 250
S   = (13.6*10**3)**2

# load pdf from lhapdf
GLUON_ID = 21
PDF = lhapdf.mkPDF("PDF4LHC15_nnlo_30_pdfas", 0)



def lumi(x1, s12, mu):
	return PDF.xfxQ2(GLUON_ID, x1, mu**2)/x1*PDF.xfxQ2(GLUON_ID, s12/(x1*S), mu**2)/s12

    

def sigma(s12, tau, alpha):
	coeff = np.pi/s12*s12**2*alpha**2/(16*v**2*(4*np.pi)**2)
	F_12 = (2/3)**2
	
	integrated_lumi = np.array(
	[integrate.quad(lumi, t, 1., args=(s, m_H))[0] for s, t in zip(s12, tau)],
	dtype=R_DTYPE
	)
	
	return coeff*F_12*integrated_lumi



def main():
	s12 = START + STEP*np.arange(NUM_POINTS, dtype=R_DTYPE)

	alpha = PDF.alphasQ2(m_H**2)
	tau   = s12/S

	sigma_ana_hT = sigma(s12, tau, alpha)

	s12.tofile('sigma_ana_hT_s12.bytes')
	sigma_ana_hT.tofile('sigma_ana_hT.bytes')



if __name__ == "__main__":
	main()
