PHIFILE='additional_files/test_data_set_lenspotentialCls.dat'
LENSEDCMBFILE='additional_files/test_data_set_lensedCls.dat'
FWHM=3.5
NOISE_LEVEL=17.7
LMIN=2
LMAXOUT=1000
LMAX=1000
LMAX_TT=1000
TMP_OUTPUT='./output'

time python LensingBiases.py -phifile ${PHIFILE} -lensedcmbfile ${LENSEDCMBFILE} -FWHM ${FWHM} -noise_level ${NOISE_LEVEL} -lmaxout ${LMAXOUT} -lmin ${LMIN} -lmax ${LMAX} -lmax_TT ${LMAX_TT} -tmp_output ${TMP_OUTPUT}
