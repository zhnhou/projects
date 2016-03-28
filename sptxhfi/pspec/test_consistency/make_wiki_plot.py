# -*- coding: utf-8 -*-
# <nbformat>4</nbformat>

# <codecell>

%matplotlib inline
%config InlineBackend.figure_format = 'svg'
import numpy as np
import matplotlib.pyplot as plt
from hpylib.util.remote_data import *
from hpylib.util.sptsz_end2end import *

# <codecell>

def do_143(end2end_files):
    dls_theory_file = '~/data_midway/projects/sptxhfi/simulations/input/dls_input_spt_150.txt'
    
    n1 = restore_end_save(sync_from_remote('midway', end2end_files[0]))
    n2 = restore_end_save(sync_from_remote('midway', end2end_files[1]))
    dls_theory = np.loadtxt(sync_from_remote('midway', dls_theory_file), usecols=[1])
    dls_theory = np.tile(dls_theory, (3,1))
    
    winfunc1 = n1['winfunc_sims']
    winfunc2 = n2['winfunc_sims']
    
    winminell = n1['winminell']
    winmaxell = n1['winmaxell']
    bands = n1['bands']

    dbs_theory1 = dls2dbs(dls_theory, winfunc1, winminell=winminell, winmaxell=winmaxell)
    dbs_theory2 = dls2dbs(dls_theory, winfunc2, winminell=winminell, winmaxell=winmaxell)
    
    tmp = n1['dbs_sims']
    dbs1_sims_ave = np.mean(tmp, axis=0)
    tmp = n2['dbs_sims']
    dbs2_sims_ave = np.mean(tmp, axis=0)
    
    plt.plot(bands, dbs1_sims_ave[1,:]/dbs_theory1[1,:], label=r'$150 \times 143$ (mask MM)')
    plt.plot(bands, dbs2_sims_ave[1,:]/dbs_theory2[1,:], label=r'$150 \times 143$ (mask 2013)')
    
    plt.plot(bands, dbs1_sims_ave[2,:]/dbs_theory1[2,:], label=r'$143 \times 143$ (mask MM)')
    plt.plot(bands, dbs2_sims_ave[2,:]/dbs_theory2[2,:], label=r'$143 \times 143$ (mask 2013)')
    
    plt.plot([0,3000],[1.0,1.0], '--', color='black')
    plt.xlim([500,3000])
    plt.ylim([0.90,1.1])
    plt.xlabel(r'$\ell$',fontsize=18)
    plt.ylabel(r'$D_{b,\rm{sims}} / D_{b, \rm{theory}}$',fontsize=18)
    
    plt.legend(loc=3)
    plt.savefig('test_150x143_143x143.pdf', format='pdf')
    plt.clf()
    plt.cla()

# <codecell>

def do_217(end2end_file):
    dls_150x220_file     = '~/data_midway/projects/sptxhfi/simulations/input/dls_ave_150x220.txt'
    dls_220_file     = '~/data_midway/projects/sptxhfi/simulations/input/dls_input_spt_220.txt'
    
    dls1_theory = np.loadtxt(sync_from_remote('midway', dls_150x220_file), usecols=[1])
    dls2_theory = np.loadtxt(sync_from_remote('midway', dls_220_file), usecols=[1])
    dls_theory = np.tile(dls1_theory, (3,1))
    dls_theory[2,:] = dls2_theory
    
    n = restore_end_save(sync_from_remote('midway', end2end_file))
    winfunc = n['winfunc_sims']
    
    winminell = n['winminell']
    winmaxell = n['winmaxell']
    bands = n['bands']
    
    dbs_theory = dls2dbs(dls_theory, winfunc, winminell=winminell, winmaxell=winmaxell)
    
    tmp = n['dbs_sims']
    dbs_sims_ave = np.mean(tmp, axis=0)
    
    plt.plot(bands, dbs_sims_ave[1,:]/dbs_theory[1,:], label=r'$150 \times 217$ (mask MM)')
    plt.plot(bands, dbs_sims_ave[2,:]/dbs_theory[2,:], label=r'$217 \times 217$ (mask MM)')
    
    plt.plot([0,3000],[1.0,1.0], '--', color='black')
    plt.xlim([500,3000])
    plt.ylim([0.90,1.1])
    plt.xlabel(r'$\ell$',fontsize=18)
    plt.ylabel(r'$D_{b,\rm{sims}} / D_{b, \rm{theory}}$',fontsize=18)
    
    plt.legend(loc=3)
    plt.savefig('test_150x217_217x217.pdf', format='pdf')
    plt.clf()
    plt.cla()

# <codecell>

s150s_h143s_files = ['~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_hfi_sn/end_combined_spt150s_hfi143s.sav',
                     '~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_hfi_sn_run_06p1/end_combined_spt150s_hfi143s.sav']
do_143(s150s_h143s_files)

# <codecell>

s150s_h217s_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_hfi_sn/end_combined_spt150s_hfi217s.sav'

do_217(s150s_h217s_file)

# <codecell>

pwd

# <codecell>


