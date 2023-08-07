import matplotlib
#matplotlib.use("TkAgg")
#import matplotlib.backends.backend_gtkagg
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
import params_class as own_classes
#import functions as func
import compute_tf as tf

# define basic class of parameters
P = own_classes.params()



number_of_time_points = 10000#200000
sample_and_integration_time_step = 0.0005#0.0001#0.0005

# initialize signal object to treat signals of dynamics
signal = own_classes.signal(  
    number_of_time_points, 
    sample_and_integration_time_step
    )
    
# initialize global dynamics object to investigate spatial mean dynamics
gd = own_classes.global_dynamics(signal)
#gd.simulation()
#quit()

#signal.compute_timefrequencyplots()
#signal.plot_tfresults()

### initialize field object to simulate network dynamics
field = own_classes.field(signal)
field.initialize_kernels()
### simulate field
field.simulation()

# plot results
#signal.plot_spacetimeresults()
#gd.plot_globalresults()


