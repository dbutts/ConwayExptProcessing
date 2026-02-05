%%
clear all
config_dir = '/mnt/isilon/code/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/Kilosort_config/';
subject = 'Sprout';
%%
Nchannels = 96;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = [2:9, 1:10, 1:10, 1:10, 1:10, 1:10, 1:8, 1, 10, 1:10, 1:10, 2:9];
xcoords   = xcoords(:);
ycoords   = [ones(1,8), ones(1,10)*2, ones(1,10)*3, ones(1,10)*4, ones(1,10)*5, ones(1,10)*6,  ones(1,8)*7, 10, 10, ones(1,10)*8, ones(1,10)*9, ones(1,8)*10];
ycoords   = ycoords(:);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 40000; % sampling frequency

fname = 'V_UT1_chanMap.mat'
save(fullfile(config_dir,subject))
save('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/V_UT1_chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%%
chanmap_Utah=nan(10,10);
for cc=1:96;
chanmap_Utah(xcoords(cc),ycoords(cc))=cc;
end
chanmap_Utah=chanmap_Utah';
%% ungrouped
Nchannels = 96;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = [2:9, 1:10, 1:10, 1:10, 1:10, 1:10, 1:8, 1, 10, 1:10, 1:10, 2:9];
xcoords   = xcoords(:);
ycoords   = [ones(1,8), ones(1,10)*2, ones(1,10)*3, ones(1,10)*4, ones(1,10)*5, ones(1,10)*6,  ones(1,8)*7, 10, 10, ones(1,10)*8, ones(1,10)*9, ones(1,8)*10];
ycoords   = ycoords(:);
kcoords   = 1:Nchannels; % grouping of channels (i.e. tetrode groups)

fs = 40000; % sampling frequency
save('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/V_UT1_chanMap_nogroup.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%
clear all
load('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/V_UT1_chanMap.mat')
cur_indx=[1 48];

connected = connected(cur_indx(1):cur_indx(2));
chanMap   = chanMap(cur_indx(1):cur_indx(2))-(cur_indx(1)-1);
chanMap0ind = chanMap0ind(cur_indx(1):cur_indx(2))-(cur_indx(1)-1);

xcoords   = xcoords(cur_indx(1):cur_indx(2));
ycoords   = ycoords(cur_indx(1):cur_indx(2));
kcoords   = kcoords(cur_indx(1):cur_indx(2));

fs = 40000; % sampling frequency
save(['/home/felix/Dropbox/Project_BevilColor/Kilosort_config/V_UT1_chanMap_' num2str(cur_indx(1)) 'to' num2str(cur_indx(2)) '.mat'], ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%
clear all
load('/home/felix/Dropbox/Project_BevilColor/Kilosort_config/V_UT1_chanMap_nogroup.mat')
cur_indx=[81 96];

connected = connected(cur_indx(1):cur_indx(2));
chanMap   = chanMap(cur_indx(1):cur_indx(2))-(cur_indx(1)-1);
chanMap0ind = chanMap0ind(cur_indx(1):cur_indx(2))-(cur_indx(1)-1);

xcoords   = xcoords(cur_indx(1):cur_indx(2))*100;
ycoords   = ycoords(cur_indx(1):cur_indx(2))*100;
kcoords   = kcoords(cur_indx(1):cur_indx(2));

fs = 40000; % sampling frequency
save(['/home/felix/Dropbox/Project_BevilColor/Kilosort_config/V_UT1_chanMap_' num2str(cur_indx(1)) 'to' num2str(cur_indx(2)) 'nogroup.mat'], ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%
