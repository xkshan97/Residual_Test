import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 
os.environ['OMP_NUM_THREADS'] = '1'
import bilby
from gwpy.timeseries import TimeSeries
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import numpy as np


sample_rate = 4096 #time series sample rate


"""This is for the Unlensed event"""
Injected_parameter_unlensed = np.loadtxt('./InjectedParameters4UnLensEvent.csv', delimiter=',')
#Reading parameters for 200 Unlensed GW events. 
#Columns from 0~11 stand for:
#Series start time, GW event index, GW duration when start frequency is 20 Hz, GW redshift, m1, m2, spin1z, spin2z, inclination, polarization, ra, dec

"""Read data"""
for event_i in range(len(Injected_parameter_unlensed)):
    start_time = int(Injected_parameter_unlensed[event_i][0])
    GWEventIndex = int(Injected_parameter_unlensed[event_i][1])
    fname_list = ['./Unlensed/H1/H-H1_GWOSC_4KHZ_R1-'+str(start_time)+'-4096_' + str(GWEventIndex) + '.gwf' , './Unlensed/L1/L-L1_GWOSC_4KHZ_R1-'+str(start_time)+'-4096_' + str(GWEventIndex) + '.gwf', './Unlensed/V1/V-V1_GWOSC_4KHZ_R1-'+str(start_time)+'-4096_' + str(GWEventIndex) + '.gwf']
    fchannel_list = ['H1:GWOSC_4KHZ_R1_STRAIN', 'L1:GWOSC_4KHZ_R1_STRAIN', 'V1:GWOSC_4KHZ_R1_STRAIN']
    psd_list = ['./PSD/H1/H1PSD.gwf' , './PSD/L1/L1PSD.gwf','./PSD/V1/V1PSD.gwf']
    psd_channel = ['H1', 'L1', 'V1']
    
    print(' duration 20 Hz = ' + str(Injected_parameter_unlensed[event_i][2]))
        
    for det_index , det in enumerate(["H1", "L1", "V1"]):#(["H1", "L1", "V1", "K1"]):
        ifo = bilby.gw.detector.get_empty_interferometer(det)
        data = TimeSeries.read(fname_list[det_index], channel = fchannel_list[det_index])
        data_start_20_index = len(data)//2
        data_end_20_index = int(data_start_20_index + sample_rate * Injected_parameter_unlensed[event_i][2])
            
        data = data[data_start_20_index: data_end_20_index] # extract the fragment where the signal exists.
        ifo.strain_data.set_from_gwpy_timeseries(data)





"""This is for the lensed event"""
Injected_parameter = np.loadtxt('./InjectedParameters4LensEvent.csv', delimiter=',')
#Reading parameters for 338 GW events. 
#Columns from 0~12 stand for:
#Series start time, GW event index, strong lensing signal order (0 or 1, where 0 means the first signal and 1 stands for the second signal), GW duration when start frequency is 20 Hz, GW redshift, m1, m2, spin1z, spin2z, inclination, polarization, ra, dec

Lens_Type_set = np.loadtxt('./LensType.csv', delimiter=',', dtype='str')
#Lens type for a specific event (Minimum or Saddle)

"""Read data"""
for event_i in range(len(Injected_parameter)):
    Lens_Type = Lens_Type_set[event_i]
    start_time = int(Injected_parameter[event_i][0])
    GWEventIndex = int(Injected_parameter[event_i][1])
    LensSignalIndex = int(Injected_parameter[event_i][2])
    fname_list = ['./' +  Lens_Type + 'WithMicro/H1/H-H1_GWOSC_4KHZ_R1-'+str(start_time)+'-4096_' + str(GWEventIndex) + "_" + str(LensSignalIndex) + '.gwf' , './' + Lens_Type + 'WithMicro/L1/L-L1_GWOSC_4KHZ_R1-'+str(start_time)+'-4096_' + str(GWEventIndex) + "_" + str(LensSignalIndex) + '.gwf', './' + Lens_Type + 'WithMicro/V1/V-V1_GWOSC_4KHZ_R1-'+str(start_time)+'-4096_' + str(GWEventIndex) + "_" + str(LensSignalIndex) + '.gwf']
    fchannel_list = ['H1:GWOSC_4KHZ_R1_STRAIN', 'L1:GWOSC_4KHZ_R1_STRAIN', 'V1:GWOSC_4KHZ_R1_STRAIN']
    psd_list = ['./PSD/H1/H1PSD.gwf' , './PSD/L1/L1PSD.gwf','./PSD/V1/V1PSD.gwf']
    psd_channel = ['H1', 'L1', 'V1']
    
    print(' duration 20 Hz = ' + str(Injected_parameter[event_i][3]))
        
    for det_index , det in enumerate(["H1", "L1", "V1"]):#(["H1", "L1", "V1", "K1"]):
        ifo = bilby.gw.detector.get_empty_interferometer(det)
        data = TimeSeries.read(fname_list[det_index], channel = fchannel_list[det_index])
        data_start_20_index = len(data)//2
        data_end_20_index = int(data_start_20_index + sample_rate * Injected_parameter[event_i][3])
            
        data = data[data_start_20_index: data_end_20_index] # extract the fragment where the signal exists.
        ifo.strain_data.set_from_gwpy_timeseries(data)
