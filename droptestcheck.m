function spk_offset = droptestcheck(plexon_fname)
disp('Drop test check starting')

pl2 = PL2ReadFileIndex(plexon_fname);

numDigitsInLastSpkChan = ceil(log10(length(pl2.SpikeChannels)));

numAIchans = sum(cellfun(@(x) contains(x.Name, 'AI'), pl2.AnalogChannels));
numDigitsInLastAIchan = ceil(log10(numAIchans));


[fs, n, ts, fn, ~] = plx_ad_v(plexon_fname,...
    ['SPKC' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);

[fs_aux, n_aux, ts_aux, fn_aux, ~] = plx_ad_v(plexon_fname, ['AI' num2str(1, ['%0' num2str(numDigitsInLastAIchan), '.f'])] );

spikeChannel1SignalDurSec = n/fs; % samples / (samples/sec)
dpiSyncSignalDurSec = n_aux/fs_aux;

droptestcheck = dpiSyncSignalDurSec - spikeChannel1SignalDurSec;
disp(['Plexon-Kofiko offset in seconds: ' num2str(droptestcheck)]) % this will tell us if the plexon time alignment issue is present

if abs(droptestcheck)>0.1
    warning("Danger - Plexon might have dropped frames! Check pl2 file.")
    spk_offset = droptestcheck;
else
    spk_offset = 0;
end

disp('Drop test check complete')


end