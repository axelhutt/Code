#!/usr/bin/env python3

from idtxl.active_information_storage import ActiveInformationStorage
#from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
from idtxl.visualise_graph import plot_network

import hde_fast_embedding

data = Data()

data.generate_mute_data(100, 5)

settings = {
    
    'cmi_estimator': 'JidtKraskovCMI',
    
    'n_perm_max_stat': 200,
    
    'n_perm_min_stat': 200,
    
    'max_lag': 5,
    
    'tau': 1
    
    }

processes = [1, 2, 3]

network_analysis = ActiveInformationStorage()
print(".....analyse_network")
results = network_analysis.analyse_network(settings, data,
    
                                           processes)