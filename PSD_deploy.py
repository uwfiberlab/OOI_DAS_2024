import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import gc  # garbage collector
import time
import h5py  # DAS usualy comes in HDF5 format
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

import scipy.signal as sgn  # for signal processing
from scipy.signal import decimate
from scipy.signal import filtfilt
from obspy import UTCDateTime  # for time conversion
from datetime import datetime  # for time conversion
from tqdm import tqdm  # progress bar

import matplotlib  # for plotting
# matplotlib.use('Agg')  # faster backend
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']  
matplotlib.rcParams['font.size'] = 20

from functools import partial  # for parallel processing
from multiprocessing import Pool, get_context  # for parallel processing


def extract_metadata(h5file, machine_name='optodas'):
    """Extract metadata from DAS HDF
    Args:
        h5file (str): path to DAS HDF file
        machine_name (str): name of interrogator
    Returns:    
        gl (float): gauge length in meters
        t0 (float): start time in seconds since 1 Jan 1970
        dt (float): sample interval in seconds
        fs (float): sampling rate in Hz
        dx (float): channel interval in meters
        un (str): unit of measurement
        ns (int): number of samples
        nx (int): number of channels
    """
    if machine_name == 'optodas':
        with h5py.File(h5file, 'r') as fp:
            gl = fp['header/gaugeLength'][()]
            t0 = fp['header/time'][()]
            dt = fp['header/dt'][()]
            fs = 1./dt
            dx = fp['header/dx'][()] # not sure why this is incorrect
            un = fp['header/unit'][()]
            ns = fp['/header/dimensionRanges/dimension0/size'][()]
            nx = fp['/header/dimensionRanges/dimension1/size'][()][0]
    elif machine_name == 'onyx':
        with h5py.File(h5file,'r') as fp:      
            gl = fp['Acquisition'].attrs['GaugeLength']
            t0 = fp['Acquisition']['Raw[0]']['RawDataTime'][0]/1e6
            fs = fp['Acquisition']['Raw[0]'].attrs['OutputDataRate']
            dt = 1./fs
            dx = fp['Acquisition'].attrs['SpatialSamplingInterval']
            un = fp['Acquisition']['Raw[0]'].attrs['RawDataUnit']
            ns  = len(fp['Acquisition']['Raw[0]']['RawDataTime'][:])
            nx = fp['Acquisition']['Raw[0]'].attrs['NumberOfLoci']
    else:
        raise ValueError('Machine name not recognized')
            

    return gl, t0, dt, fs, dx, un, ns, nx


### functions to calculate PDF of multiple channels
### Modified from Enthan Williams's code
def ppsd(data,fs,fmin,fmax):
    """
    data:  2D array, the statistics is calculated along axis=0
    fs: sampling rate
    fmin: minimum frequency for statistics
    fmax: maximum frequency for statictics
    """
    ns = data.shape[1]
    nx = data.shape[0]
    
    ### Demean, detrend
    data -= np.mean(data, axis=1, keepdims=True) 
    data = sgn.detrend(data, axis=1) 
    
    freq, spec = sgn.periodogram(data, fs, window='hamming', axis=-1)

    # spec = np.abs(np.fft.fft(data))
    # freq = np.fft.fftfreq(data.shape[1], 1/fs)

    freq = np.tile(freq,(nx,1)).flatten()

    # print(freq.shape, spec.shape)


    ### Generate PDF
    xbins = np.logspace(np.log10(fmin),np.log10(fmax),100)
    ybins = np.logspace(np.log10(np.nanmin(spec)),np.log10(np.nanmax(spec)),200)
    
    H,xe,ye = np.histogram2d(freq.flatten(), spec.flatten(), bins=(xbins,ybins))
  
    
    return H/np.nansum(H, axis=1, keepdims=True), (xe[1:] + xe[:-1])/2, (ye[1:] + ye[:-1])/2
    

def psd_stats(H,xm,ym):   
    ym = np.log10(ym)
    mean = np.zeros(len(xm))
    variance = mean.copy()
    for ix in range(len(xm)):
        mean[ix] = np.average(ym,weights=H[ix,:])
        variance[ix] = np.average((ym-mean[ix])**2,weights=H[ix,:])
    
    return xm,10**mean,variance


def read_decimate(file_path, dsamp_factor=20, start_ch=0, end_ch=100, machine_name='onyx'):
    if machine_name == 'optodas':
        with h5py.File(file_path, 'r') as f:
            minute_data = f['data'][:, start_ch:end_ch].T
    elif machine_name == 'onyx':
        with h5py.File(file_path,'r') as f:      
            minute_data = f['Acquisition']['Raw[0]']['RawData'][:, start_ch:end_ch].T
    else:
        raise ValueError('Machine name not recognized')
    
    if dsamp_factor>1:
        downsample_data = decimate(minute_data, q=dsamp_factor, ftype='fir', zero_phase=True)   
    else:
        downsample_data = minute_data
    
    return downsample_data


def ppsd_on_fly(data_dir,out_dir, machine_name, num_proc, start_ch, end_ch, start_file, num_file, 
dsamp_factor, dx_correct,channel_bin, channel_interval, correction_factor, amp_type='strain_rate', on_macos=False):
    
    file_list = []
    with os.scandir(data_dir) as entries:
        for entry in sorted(entries, key=lambda e: e.name)[start_file:start_file + num_file]:
            file_list.append(entry.name)

    print('#'*10 + ' Working on '+str(file_list[0])+' to '+str(file_list[-1])+ ' ' + '#'*10)

    file_path = [os.path.join(data_dir,i) for i in file_list]

    gl, t0, dt, fs, dx, un, ns, nx = extract_metadata(file_path[0], machine_name='optodas')

    dx = dx * dx_correct

    new_fs = int(fs / dsamp_factor)        # final sample rate after downsampling

    # %% multi-process to read and decimate lots of files 
    partial_func = partial(read_decimate, dsamp_factor=dsamp_factor, start_ch=start_ch, end_ch=end_ch, machine_name=machine_name)
    if on_macos:
        with get_context('fork').Pool(processes=num_proc) as pool:   # pool is closed automatically and join as a list
            print("# threads, using fork context on Mac: ", num_proc)
            full_time = pool.map(partial_func, file_path[:])
    else:
        with Pool(processes=num_proc) as pool:
            print("# threads: ", num_proc)
            full_time = pool.map(partial_func, file_path[:])

    # %% concatenate the list elements in time
    print("- Concatenating time series")
    full_time_data = np.concatenate(full_time, axis=1)

    print(f'final shape: {full_time_data.shape}')

    ### Convert phase to strain/ strain rate
    hour_data = full_time_data * correction_factor

    ### Plot PSD
    print("- Ploting PSD")
    fig,ax = plt.subplots(1,1,figsize=(20,10))

    ### divide into several segments, indicated by different colors
    half_bin = int(channel_bin/2)
    ch_list = np.arange(half_bin+1,(end_ch-start_ch),channel_interval)[::-1]
    colors = matplotlib.colormaps['Blues'](np.linspace(0.15, 1, len(ch_list)))[::-1]

    

    for i, chan in tqdm(enumerate(ch_list)):
        
        trs = hour_data[chan-half_bin:chan+half_bin,:]
        if amp_type == 'strain':
            trs = np.diff(trs, axis=1) * new_fs

        H,xm,ym = ppsd(trs,new_fs,2e-4,1e2)
        xm,mn,vr = psd_stats(H,xm,ym)

        ax.plot(1/xm, mn,linewidth=5, label='%.1f m' % (chan*dx), zorder=2, color=colors[i])
        if i == 0:
            img=ax.pcolormesh(1/xm,ym,H.T,cmap='hot_r',vmin=0,vmax=0.1, zorder=1)

    ax.set_xscale('log')
    ax.set_yscale('log')   
    ax.set_ylim([1e-20,1e-9])
    ax.set_xlabel('Period (s)')
    ax.set_ylabel('PSD rel. strain rate ^2')
    ax.grid(which='major', color='#DDDDDD', linewidth=3, zorder=0)
    ax.grid(which='minor', color='#EEEEEE', linewidth=2, linestyle='--', zorder=0)
    plt.colorbar(img, ax=ax, aspect=50).set_label('Probability')
    plt.legend(loc="upper left")
    plt.tight_layout()  
    plt.savefig(out_dir + str(file_list[0])+str(file_list[-1])+'PSD.png',dpi=300)


if __name__ == '__main__':

    for starting_id in np.arange(0, 1, 1):
        since = time.time()
        ppsd_on_fly(data_dir = '/Users/qb/Downloads/Archive/data_h5', 
                    out_dir = '/Users/qb/Downloads/Archive/',
                    machine_name = 'optodas', 
                    num_proc = 2, 
                    start_ch = 0, 
                    end_ch = 9000, 
                    start_file = starting_id, 
                    num_file = 20, 
                    dsamp_factor = 1, 
                    dx_correct = 5,
                    channel_bin = 200, 
                    channel_interval = 2000, 
                    correction_factor = (1550.12 * 1e-9) / (0.78 * 4 * np.pi * 1.4677), 
                    amp_type='strain_rate',
                    on_macos=True)

        time_elapsed = time.time() - since
        print('- Time taken: {:.0f}m {:.0f}s'.format(time_elapsed // 60, time_elapsed % 60))
