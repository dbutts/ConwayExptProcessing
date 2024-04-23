%%
%%

Nchannels = 28;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = [1 1 1 2 1 2 3 1 1 3 1 4 1 1 2 1 2 3 4 5 6 1 2 1 1 2 1 1 ];
ycoords   = [1 2 3 3 4 4 4 5 6 6 7 8 9 10 10 11 11 11 11 11 11 13 13 12 16 16 15 14];
xcoords   = xcoords(:);
ycoords   = ycoords(:);
kcoords   = 1:Nchannels; % grouping of channels (i.e. tetrode groups)
%%
fs = 40000; % sampling frequency
save('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Nform_chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
