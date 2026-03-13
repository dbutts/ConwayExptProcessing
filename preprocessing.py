import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.exporters as sexp
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from kilosort import io
import pandas as pd
import neo.rawio
import os

def makeProbeFromFile(chanMapDict, grouped=0, chnOffset =0):
    chanMap = pd.read_excel(chanMapDict['filename'],sheet_name = chanMapDict['sheet']).iloc[chanMapDict['cellIdx'][0], chanMapDict['cellIdx'][1]]
    chanMap = np.rot90(np.flipud(chanMap),-1)
    [xc, yc] = [np.where(pd.notna(chanMap))[1], np.where(pd.notna(chanMap))[0]]
    
    n_chan = len(xc)
    connected = np.ones((n_chan,1)).astype('bool')
    
    chanMap = np.arange(n_chan) + chnOffset
    
    if grouped:
        kcoords = np.ones(n_chan)
    else:
        kcoords = np.arange(n_chan)

    probe = {
        'chanMap': chanMap,
        'xc': xc,
        'yc': yc,
        'kcoords': kcoords,
        'n_chan': n_chan
    }
    return probe

def makeProbeFromParameters(n_chan, spacing, grouped=0, chnOffset= 0, x0=0, y0=0):
    chanMap = np.arange(n_chan) + chnOffset
    xc = np.array([x0 + spacing[0]*i for i in reversed(range(n_chan))])
    yc = np.array([y0 + spacing[1]*i for i in reversed(range(n_chan))]) # reversed makes things identical to how they were in matlab
    
    if grouped:
        kcoords = np.ones(n_chan)
    else:
        kcoords = np.arange(n_chan)
    
    
    probe = {
        'chanMap': chanMap,
        'xc': xc,
        'yc': yc,
        'kcoords': kcoords,
        'n_chan': n_chan
    }
    return probe
    

monkey_name = "Sprout"
exptDate = "260205"

dirpath = Path.home() / "Data" / "V1_Fovea" / monkey_name/ exptDate
processingPath = Path.home() / "Git" / "ConwayExptProcessing"

# find plexon filename

plexon_fname = dirpath / (exptDate + "_153540_" + monkey_name + ".pl2")

assert plexon_fname.is_file()

#fname = "/home/bizon/Data/V1_Fovea/Sprout/260205/260205_153540_Sprout.pl2"

chanMapPath = processingPath / "chanMapFiles" / monkey_name


chanMapDicts = [{'filename': chanMapPath / 'SN 7623-000080 Mapping and Impedance' / 'SN 7623-000080.xlsm',
                'sheet': 'Z check',
                'cellIdx': [np.arange(14, 25), np.arange(15,25)]},
               {'filename': chanMapPath / 'SN 7624-000191 Mapping and Impedance' /'1125-15 SN 7624-000191.xlsm',
               'sheet': 'Z check',
               'cellIdx':[np.arange(14, 25), np.arange(15,25)]}, {}]



arrayDict = {'labels': ['UT1', 'UT2', 'lam'],
            'n_chan': [(), (), 64],
            'spacing': [(), (), (0,50)],
             'grouped': [0,0,1],
            'chanMapDict': chanMapDicts}

# for each channel, determine if it has a channel map and if not specifiy geometry
probeDicts = []
chnOffset = 0

for i in range(len(arrayDict['labels'])):
    # if this array has an excel file, make probe info from file
    if arrayDict['chanMapDict'][i]:
        probeDicts.append(makeProbeFromFile(arrayDict['chanMapDict'][i], arrayDict['grouped'][i], chnOffset)) 
        chnOffset = chnOffset + probeDicts[i]['n_chan']
    else:
        probeDicts.append(makeProbeFromParameters(arrayDict['n_chan'][i], arrayDict['spacing'][i], arrayDict['grouped'][i], chnOffset))
        chnOffset = chnOffset + probeDicts[i]['n_chan']
        
        
        
rec_SPKC = se.read_plexon2(plexon_fname, stream_id = 'SPKC')

lastSPKC = len(rec_SPKC.channel_ids)
numDigitsInLastSPKC = int(np.ceil(np.log10(lastSPKC)))

channel_ids_UT1 = ["SPKC" + f"{(n+1):0{numDigitsInLastSPKC}d}" for n in probeDicts[0]['chanMap']]

rec_SPKC_UT1 = rec_SPKC.select_channels(channel_ids= channel_ids_UT1)

job_kwargs = {
    'n_jobs': 1,          # Use 4 CPU cores
    'chunk_duration': '10s', # Process 10 second chunks at a time
    'progress_bar': True  # Display a progress bar
}


n_channels = rec_SPKC_UT1.get_num_channels()
n_frames = rec_SPKC_UT1.get_num_frames()
n_frames = 1000
dtype = 'int16'

output_file = dirpath / '260205_153540_Sprout' / 'kilosort' / 'kilosort_UT1.bin'
os.makedirs(os.path.dirname(output_file), exist_ok=True)

fp = np.memmap(output_file, dtype=dtype, mode = 'w+', shape = (n_channels, n_frames))
chunk_size = 1_000_000

for start in range(0, n_frames, chunk_size):
    end = min(start+chunk_size, n_frames)
    traces = rec_SPKC_UT1.get_traces(start_frame=start, end_frame=n_frames).astype(dtype)
    #fp[:, start:end] = traces.T
    fp[:, 0:n_frames] = traces.T # for testing
    print(f"Wrote frames {start} to {end}")

fp.flush()
print("Done! Binary file is ready for Kilosort 4")