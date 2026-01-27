%%
%%

Nchannels = 96;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = [2:9, 1:10, 1:9, 10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 2,3,1, 5:9];
xcoords   = xcoords(:);
ycoords   = [ones(1,8), ones(1,10)*2, ones(1,9)*3, 1, ones(1,10)*4, ones(1,10)*5, ones(1,10)*6, ones(1,10)*7, ones(1,10)*8, ones(1,10)*9, ones(1,8)*10];
ycoords   = ycoords(:);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 40000; % sampling frequency
save('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/Utah_chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')