%%
clear all
subject = 'Sprout';
grouped = 0;

config_dir = '/mnt/isilon/code/ConwayExptProcessing/Dependencies/Kilotools_FB_2023/Kilosort_config/';


%7623

%Electrode number viewing from pad side
chanMap_Utah = [nan	88	78	68	58	48	38	28	18	nan
                96	87	77	67	57	47	37	27	17	8
                95	86	76	66	56	46	36	26	16	7
                94	85	75	65	55	45	35	25	15	6
                93	84	74	64	54	44	34	24	14	5
                92	83	73	63	53	43	33	23	13	4
                91	82	72	62	52	42	32	22	12	3
                90	81	71	61	51	41	31	21	11	2
                89	80	70	60	50	40	30	20	10	1
                49	79	69	59	nan	39	29	19	9	nan];

chanMap_Utah = rot90(flipud(chanMap_Utah));
[xcoords,ycoords] = find(~isnan(chanMap_Utah));


%%
Nchannels = sum(~isnan(chanMap_Utah), 'all');
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

if grouped
    kcoords = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
    fname = 'V1_UT1_chanMap.mat';
else
    kcoords = 1:Nchannels; % grouping of channels (i.e. tetrode groups)
    fname = 'V1_UT1_chanMap_nogroup.mat';
end

fs = 40000; % sampling frequency
save(fullfile(config_dir,subject, fname), 'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%%
