function stas = get_sta( data, targets, use_inds, to_shift, stim_deltas, num_lags )
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
    to_shift = 1; 
end

if nargin < 5
    NT=size(data.ETtrace,2);
    stim_deltas = zeros(2,NT); 
end

if nargin < 6
    num_lags = 6;
end

stim_shift=permute(data.stim,[4 1 2 3]);

if to_shift~=0
    stim_shift = shift_stim( stim_shift, data.ETtrace, stim_deltas );
    shiftstr = 'with shifts';
else
    shiftstr = 'no shifts';
end

%reshaping stimuli to get STAs with matrix multiplication, which is faster
stim2 = single(reshape(stim_shift,size(stim_shift,1),3*60*60))./127;

stas=zeros(max(targets),num_lags,3*60*60);

disp('Now plotting!')
for cc=targets
	tic
	figure;
	for curlag=1:num_lags
		cur_STA1(curlag,:) = binned_SU(use_inds+curlag,cc)' * stim2(use_inds,:);
	end

    stas(cc,curlag, :)=cur_STA1;
	for curlag=1:num_lags
		cur_StA2=reshape(cur_STA1(curlag,:),60,180);
		%curlim=150;%
		curlim=max(abs(cur_STA1(:)'));
		cur_StA2(:,[60 120])=curlim;

		subplot(6,1,curlag)
		imagesc(cur_StA2); clim([-curlim curlim]); pbaspect([3 1 1])

		colormap(gray); xlabel(['Lag ' num2str(curlag)]);
		if curlag == 1
			ylabel('S          L-M          Lum'); 
            if cc<=nSU;
                titlestr = sprintf('SU %0d: ch# %d; %s', cc, spk_ch(cc), shiftstr);
            else
                titlestr = sprintf('MU %0d: ch# %d; %s', cc-nSU, spk_ch(cc), shiftstr);
            end
			%title(['SU ' num2str(cc) ': ch#', num2str(spk_ch(cc)), 'spkID' num2str(spk_ID(cc))]); 
            title(titlestr)
		end
	end
	%    figtitle(['Probe ' num2str(Robs_probe_ID(cc)) '   Unit ' num2str(cc) ])
	toc
	%    pause
end

