function stas = get_sta( data, targets, use_inds, to_shift, stim_deltas, num_lags, save_vars )
%
% Usage: stas = get_sta( data, <targets>, <to_shift>, <stim_deltas> )
%
% PLAN COMING BACK:
% -- Make STA function that uses the shifted stim and spikes (and valid data)
%    It then calls the STA function, as well as shifting function, and reshaping and so on 

%% Now we process the data to make STAs
binned_SU = [single(data.Robs'), single(data.RobsMU')];
spk_ch = [data.Robs_probe_ID; data.RobsMU_probe_ID];
nSU=size(data.Robs,1); 


if nargin < 2 
    targets = 1:size(binned_SU,2); % default to plotting all SUs and MUs
end

if nargin <3 || isempty(use_inds)
    use_inds = data.valid_data;
    use_inds(end-10:end) = []; %cut last few indices to avoid artifacts
end

if nargin < 4 || isempty(to_shift)
    to_shift = 0; 
end

if nargin < 5
    NT=size(data.ETtrace,2);
    stim_deltas = zeros(2,NT); 
end

if nargin < 6 || isempty(num_lags)
    num_lags = 6;
end

if nargin < 7 || isempty(save_vars)
    save_vars.to_save=0;
end

stim_shift=permute(data.stim,[4 1 2 3]);

if to_shift==1
    stim_shift = shift_stim( stim_shift, data.ETtrace, stim_deltas );
    shiftstr = 'with shifts';
elseif to_shift==0
    stim_shift = shift_stim( stim_shift, zeros(size(data.ETtrace)), stim_deltas );
    shiftstr = 'no shifts';
else
    shiftstr = 'no shifts';
end

%reshaping stimuli to get STAs with matrix multiplication, which is faster
stim2 = single(reshape(stim_shift,size(stim_shift,1),3*60*60))./127;

stas=zeros(max(targets),num_lags,3*60*60);

disp('Now plotting!')
for cc=targets
	tic
    % cameron added 2 lines: 
    indexed_spikes = binned_SU(use_inds,cc);
    total_spikes = sum(indexed_spikes == 1);

	STA = figure;
	for curlag=1:num_lags
		cur_STA1(curlag,:) = binned_SU(use_inds+curlag,cc)' * stim2(use_inds,:);
	end
	curlim=max(abs(cur_STA1(:)')); 
        if curlim==0; curlim=0.1; end % avoids plotting bugs if a bad STA is included in a large set of plots

    stas(cc,:, :)=cur_STA1;
	for curlag=1:num_lags
		cur_StA2=reshape(cur_STA1(curlag,:),60,180);
		%curlim=150;%
		cur_StA2(:,[60 120])=curlim; % plots vertical lines dividing up the plots

		subplot(6,3,1+(curlag-1)*3)
		imagesc(cur_StA2); clim([-curlim curlim]); pbaspect([3 1 1])

        cur_StA3=reshape(cur_STA1(curlag,:),60*60,3);
        cur_StA_L = cur_StA3*[-0.040, 0.5220,0]'; % from 01_2022 calib
		cur_StA_L2=reshape(cur_StA_L,60,60);
        cur_StA_M = cur_StA3*[0.0351, -0.4625, 0]'; % from 01_2022 calib
		cur_StA_M2=reshape(cur_StA_M,60,60);

		if curlag == 1
			xlabel('Lum          L-M          S'); 
            if cc<=nSU;
                titlestr = sprintf('SU %0d: ch# %d; %s', cc, spk_ch(cc), shiftstr);
            else
                titlestr = sprintf('MU %0d: ch# %d; %s', cc-nSU, spk_ch(cc), shiftstr);
            end
			%title(['SU ' num2str(cc) ': ch#', num2str(spk_ch(cc)), 'spkID' num2str(spk_ID(cc))]); 
            title(titlestr)
            % cameron added 1 line:
            subtitle("Total exp time: " + length(use_inds)/240 + " seconds" + newline + "end" ...
                + newline + "total spikes: " + total_spikes);
        end

        subplot(6,3,2+(curlag-1)*3)
		imagesc(cur_StA_L2); clim([-curlim curlim].*.522); pbaspect([1 1 1]); % color limit based on L-isolating max
		if curlag == 1; xlabel('L'); end

		subplot(6,3,3+(curlag-1)*3)
		imagesc(cur_StA_M2); clim([-curlim curlim].*.4625); pbaspect([1 1 1])
		
		colormap(gray); xlabel(['Lag ' num2str(curlag)]);

        if curlag == 1; xlabel('M'); end


	end
	%    figtitle(['Probe ' num2str(Robs_probe_ID(cc)) '   Unit ' num2str(cc) ])
	toc
	%    pause
    if save_vars.to_save==1
        % save the STA as pdf saveas(figure, filename)
    	saveas(STA,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_STA.pdf'])
        % save the STA as a jpg
        % saveas(STA,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_STA.jpg'])
    end
end

