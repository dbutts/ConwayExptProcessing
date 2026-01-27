function stas = shifted_sta( stim, binned_SU, ETtrace, valid_data, stim_deltas )
%
% Usage: stas = shifted_sta( stim, binned_SU, ETtrace, valid_data, <stim_deltas> )
%
% PLAN COMING BACK:
% -- Make STA function that uses the shifted stim and spikes (and valid data)
% -- Make an overarching function that takes the 'data' and takes what it needs
%    It then calls the STA function, as well as shifting function, and reshaping and so on 

%% Now we process the data to make STAs
binned_SU = [single(Robs'), single(RobsMU')];
stim_shift = permute(stim,[4 1 2 3]);

%to shift the STAs with the eye-trace data
ETshifts = mod(round(-ETtrace)+30,60);
ETshifts(isnan(ETshifts))=30;

tic
for i=1:size(stim_shift,1)
	stim_shift(i,:,:,:) = circshift(stim_shift(i,:,:,:), ETshifts(1,i),2);
	stim_shift(i,:,:,:) = circshift(stim_shift(i,:,:,:), ETshifts(2,i),3);
end
toc
%end shifting of STAs with eye trace

%reshaping stimuli to get STAs with matrix multiplication, which is faster
stim2 = single(reshape(stim_shift,size(stim_shift,1),3*60*60)) ./ 127;

% stimET=permute(stimET,[4 1 2 3]); %stimX=stim;
% stim2=(double(reshape(stimET,size(stimET,1),3*60*60))./127);

%cur_use_inds = valid_data(1:100000);
cur_use_inds = valid_data;
cur_use_inds(end-10:end) = [];

disp('All done with packaging!')


% % STA STUFF
for cc=1:size(binned_SU,2)
	tic
	figure;
	for curlag=1:6
		cur_STA1(curlag,:) = binned_SU(cur_use_inds+curlag,cc)' * stim2(cur_use_inds,:);
	end

	for curlag=1:6 
		cur_StA2 = reshape(cur_STA1(curlag,:),60,180);
		%curlim=150;%
		curlim=max(abs(cur_STA1(:)'));
		cur_StA2(:,[60 120]) = curlim;

		subplot(6,1,curlag)
		imagesc(cur_StA2); clim([-curlim curlim]); pbaspect([3 1 1])

		colormap(gray); xlabel(['Lag ' num2str(curlag)]);
		if curlag == 1
			ylabel('S          L-M          Lum'); 
			title(['ch#', num2str(spk_ch(cc)), 'spkID' num2str(spk_ID(cc))]); 
		end
	end
	%    figtitle(['Probe ' num2str(Robs_probe_ID(cc)) '   Unit ' num2str(cc) ])
	toc
	%    pause
end
