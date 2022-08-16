#!/usr/bin/env python3
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
m_t = 172
v   = 250
S   = (13.6*10**3)**2

# load pdf from lhapdf
GLUON_ID = 21
PDF = lhapdf.mkPDF("PDF4LHC15_nnlo_30_pdfas", 0)



def lumi(x1, s12, mu):
    return PDF.xfxQ2(GLUON_ID, x1, mu**2)/x1*PDF.xfxQ2(GLUON_ID, s12/(x1*S), mu**2)/s12



def C_0(s12):
    tau = 4.*m_t**2/s12

    case_tau_geq_1 = (tau >= 1.0)
    case_tau_l_1   = np.logical_not(case_tau_geq_1)

    tau_geq_1 = tau[case_tau_geq_1]
    tau_l_1   = tau[case_tau_l_1]

    C_0_arr                 = np.empty(NUM_POINTS, dtype=C_DTYPE)
    C_0_arr[case_tau_geq_1]   = -2./s12[case_tau_geq_1]*np.arcsin(np.sqrt(1./tau_geq_1))**2
    C_0_arr[case_tau_l_1] = 1./(2.*s12[case_tau_l_1])*(np.log((1+np.sqrt(1.-tau_l_1))/(1.-np.sqrt(1.-tau_l_1)))-1j*np.pi)**2

    return C_0_arr



def sigma(s12, tau, alpha):
    coeff = np.pi*m_t**2*alpha**2/(16*(4*np.pi)**2*v**2*s12)
    F_12  = np.abs(2*m_t*(2+(4*m_t**2-s12)*C_0(s12)))**2

    integrated_lumi = np.array(
        [integrate.quad(lumi, t, 1., args=(s, m_H))[0] for s, t in zip(s12, tau)],
        dtype=R_DTYPE
    )

    return coeff*F_12*integrated_lumi



def main():
    s12   = START + STEP*np.arange(NUM_POINTS, dtype=R_DTYPE)
    alpha = PDF.alphasQ2(m_H**2)
    tau   = s12/S

    sigma_ana = sigma(s12, tau, alpha)

    s12.tofile('sigma_ana_s12.bytes')
    sigma_ana.tofile('sigma_ana.bytes')



if __name__ == "__main__":
    main()
