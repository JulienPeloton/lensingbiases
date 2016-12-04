PHIFILE='additional_files/planck_lensing_wp_highL_bestFit_20130627_lenspotentialCls.dat'
LENSEDCMBFILE='additional_files/planck_lensing_wp_highL_bestFit_20130627_lensedCls.dat'
FWHM=7.0
NOISE_LEVEL=27.0
LMIN=2
LMAXOUT=1000
LMAX=1000
LMAX_TT=1000
TMP_OUTPUT='./output'

time python LensingBiases.py -phifile ${PHIFILE} -lensedcmbfile ${LENSEDCMBFILE} -FWHM ${FWHM} -noise_level ${NOISE_LEVEL} -lmaxout ${LMAXOUT} -lmin ${LMIN} -lmax ${LMAX} -lmax_TT ${LMAX_TT} -tmp_output ${TMP_OUTPUT}
