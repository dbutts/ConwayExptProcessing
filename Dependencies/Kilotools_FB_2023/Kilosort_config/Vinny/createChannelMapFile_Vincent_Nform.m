%%
clear all

%%

Nchannels = 128;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = [ones(1,8)*4, ones(1,8)*3, ones(1,8)*2, ones(1,8), ...
    ones(1,8)*8, ones(1,8)*7, ones(1,8)*6, ones(1,8)*5, ...
    ones(1,8)*12, ones(1,8)*11, ones(1,8)*10, ones(1,8)*9, ...
    ones(1,8)*16, ones(1,8)*15, ones(1,8)*14, ones(1,8)*13];
ycoords   = [repmat([8:-1:1],1,16)];
xcoords   = xcoords(:);
ycoords   = ycoords(:);
kcoords   = 1:Nchannels; % grouping of channels (i.e. tetrode groups)
%%
fs = 40000; % sampling frequency
%save('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/V_NF_chanMap.mat', ...
%    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
save('/home/conwaylab/Git/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/Kilosort_config/Vinny/V_NF_chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')


%%
chanmap_NF=nan(16,8);
for cc=1:128;
chanmap_NF(xcoords(cc),ycoords(cc))=cc;
end
chanmap_NF=chanmap_NF';

%%