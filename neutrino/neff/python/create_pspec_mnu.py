import os
from hpylib.cosmology.get_pspec_param import *

home_path = os.getenv('HOME')+'/'
host_name = os.getenv('ENV_HOSTNAME')

sample_path = home_path+'data_'+host_name+'/projects/neutrino/neff/param_sample/'

#sample_prefix = 'params_mnu_const_omegab_omegac_thetas'
sample_prefix = 'params_mnu_const_omegab_omegam'
param = get_params(sample_path, sample_prefix)
#output_root = sample_path+'Dl_mnu_const_omegab_zeq_thetas'
output_root = sample_path+'Dl_mnu_const_omegab_omegam'
dl_sample = calc_pspec(param, output_root=output_root, TCMB=2.7255e6)
