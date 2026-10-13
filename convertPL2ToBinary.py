# convertPL2ToBinary
# Written by Max Greene, 4/2026


# Required arguments:
# stream_id (e.g. 'SPKC')
# requested_chan_nums (one-indexed channel numbers as they appear in plexon file)
# plexon file (full path)
# output file

import openreadout
import numpy as np
from pathlib import Path
import os


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
requested_chan_nums = parseChannelRanges(input("Enter channels/channel ranges (e.g. 1:32): "))
array_name = input("Enter array label (e.g. UT1, laminar): ")

plexon_dir = os.path.dirname(plexon_fname)
plexon_basename = os.path.basename(plexon_fname)
plexon_basename_no_ext,_ = os.path.splitext(plexon_basename)

dtype = np.int16

with openreadout.open(plexon_fname) as f:
    trace_names = [trace['name'] for trace in f.traces]
    trace_idx = trace_names.index(stream_id)

    channels = f.traces[trace_idx]['channels']
    all_chan_names = [x['name'] for x in channels]

    requested_chan_names = [
        chname for chname in all_chan_names
        if int(chname.removeprefix(stream_id)) in requested_chan_nums
    ]

    requested_chan_indices = np.array([
        all_chan_names.index(chname)
        for chname in requested_chan_names
    ])

    n_samples = f.traces[trace_idx]['sample_count']
    chunk_size = 10_000_000

    dat_basename = plexon_basename_no_ext + '_' + array_name + '_' + f"{min(requested_chan_nums)}to{max(requested_chan_nums)}"

    dat_fname = os.path.join(
        plexon_dir,
        plexon_basename_no_ext,
        'kilosort' + '_' + array_name + '_' + f"{min(requested_chan_nums)}to{max(requested_chan_nums)}",
        dat_basename + '.dat'
    )

    os.makedirs(os.path.dirname(dat_fname), exist_ok=True)

    out = np.memmap(
        dat_fname,
        dtype=dtype,
        mode="w+",
        shape=(n_samples, len(requested_chan_nums))
    )

    print("Output file created")

    scale_by_name = {ch['name']: ch['scale'] for ch in channels}
    scales = np.array([
        scale_by_name[ch]
        for ch in requested_chan_names
    ])

    for first_sample in range(0, n_samples, chunk_size):
        last_sample = min(first_sample + chunk_size, n_samples)

        traces = f.read_trace(
            trace=trace_idx,
            first_sample=first_sample,
            max_samples=last_sample - first_sample
        )

        traces = np.asarray(traces)
        traces = traces[requested_chan_indices, :]

        out[first_sample:last_sample, :] = np.rint(
            traces / scales[:, None]
        ).T.astype(np.int16)

        print(f"{last_sample} / {n_samples} samples complete")

    out.flush()
    del out

print("Done! Binary file is ready for Kilosort 4")