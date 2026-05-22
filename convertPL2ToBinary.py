# convertPL2ToBinary
# Written by Max Greene, 4/2026



# Required arguments:
# stream_id (e.g. 'SPKC')
# channel_nums (one-indexed channel numbers as they appear in plexon file)
# plexon file (full path)
# output file

# parse channel number input


import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.exporters as sexp
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import neo.rawio
import os
import spikeinterface as si
import psutil


def parseChannelRanges(s):
    result = []
    for part in s.split(','):
        part = part.strip()
        if ':' in part:
            start, end = map(int, part.split(':'))
            result.extend(range(start, end+1))
        else:
            result.append(int(part))
    return result

plexon_fname = os.path.abspath(input("Enter full plexon file name: "))
assert os.path.isfile(plexon_fname)
stream_id = input("Enter channel stream id (e.g. SPKC): ")
channel_nums = parseChannelRanges(input("Enter channels/channel ranges (e.g. 1:32): "))
array_name = input("Enter array label (e.g. UT1, laminar): ")

plexon_dir = os.path.dirname(plexon_fname)
plexon_basename = os.path.basename(plexon_fname)
plexon_basename_no_ext,_ =  os.path.splitext(plexon_basename) 

# based on available RAM and disk space, figure out how much of the file can be safely read:
# make sure arrays are at most ~50-70% of free RAM, and memmap is ~80-90% of free disk space

mem = psutil.virtual_memory()
usage = psutil.disk_usage('Z:\\')

max_array_size_bytes = int(0.6*mem.available)
max_array_int16_samples = max_array_size_bytes // 2

max_memmap_bytes = int(0.85*usage.free)
max_memmap_int16_samples = max_memmap_bytes // 2

# monkey_name = "Sprout"
# exptDate = "260205"
#pl2path = Path("Z:/") / monkey_name/ exptDate
#configPath = Path("K:/") / monkey_name 
# find plexon filename
#plexon_fname = pl2path / (exptDate + "_153540_" + monkey_name + ".pl2")

print("Getting recording info")
recording = se.read_plexon2(plexon_fname, stream_id = stream_id)

fs = recording.get_sampling_frequency()
n_frames = recording.get_num_frames()
lastSPKC = len(recording.channel_ids)
numDigitsInLastSPKC = int(np.ceil(np.log10(lastSPKC)))
channel_ids = ["SPKC" + f"{(n):0{numDigitsInLastSPKC}d}" for n in channel_nums]

chunk_size = np.min([max_array_int16_samples, n_frames])

dtype = np.int16

max_memmap_shape = (n_frames, max_memmap_int16_samples // n_frames)

n_channels = np.min([len(channel_nums), max_memmap_shape[1]])

if n_channels != len(channel_nums):
    print(f"WARNING: Processing {len(channel_nums)} channels may cause memory error. {n_channels} will be processed.")

dat_basename =  plexon_basename_no_ext + '_' + array_name + '_' + f"{min(channel_nums)}to{min(channel_nums) + n_channels-1}"
dat_fname = os.path.join(plexon_dir, plexon_basename_no_ext, 'kilosort' + '_' + array_name + '_' + f"{min(channel_nums)}to{min(channel_nums) + n_channels-1}", dat_basename + '.dat')
os.makedirs(os.path.dirname(dat_fname), exist_ok=True)
out = np.memmap(dat_fname, dtype=dtype, mode="w+", shape=(n_frames, n_channels))

print("Output file created")

start_chan = 0
for ch in channel_nums:
    chan_id = channel_ids[ch-1]

    for start_frame in range(0, n_frames, chunk_size):
        end_frame = min(start_frame + chunk_size, n_frames)

        traces = recording.get_traces(
            channel_ids=[chan_id],
            start_frame=start_frame,
            end_frame=end_frame,
        )

        out[start_frame:end_frame, ch-1] = traces[:, 0].astype(dtype)
        del traces

    print(f"{chan_id} read complete")

out.flush()
del out
        
print("Done! Binary file is ready for Kilosort 4")
