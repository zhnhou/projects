import os
from hpylib.cosmology.get_pspec_param import *

home_path = os.getenv('HOME')+'/'
host_name = os.getenv('ENV_HOSTNAME')

sample_path = home_path+'data_'+host_name+'/projects/neutrino/neff/param_sample/'

sample_prefix = 'params_base_TT_lowP_lensing_const_omegab_zeq_thetas'
param = get_params(sample_path, sample_prefix)
output_root = sample_path+'Dl_base_TT_lowP_lensing_const_omegab_zeq_thetas'
dl_sample = calc_pspec(param, output_root=output_root)

sample_prefix = 'params_base_TT_lowP_lensing_const_omegab_zeq_thetas_thetad'
param = get_params(sample_path, sample_prefix)
output_root = sample_path+'Dl_base_TT_lowP_lensing_const_omegab_zeq_thetas_thetad'
dl_sample = calc_pspec(param, output_root=output_root)
