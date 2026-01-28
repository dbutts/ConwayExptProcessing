function [fwdc_L, fwdc_M, fwdc_diff] = get_fwdc( data, targets, use_inds, to_shift, stim_deltas, thresh, num_lags, save_vars )
%
% Usage: fwdc = get_fwdc( data, <targets>, <to_shift>, <stim_deltas> )
%

%% First we process the data just like when we make STAs
binned_SU = [single(data.Robs'), single(data.RobsMU')];
spk_ch = [data.Robs_probe_ID; data.RobsMU_probe_ID];
nSU=size(data.Robs,1); 
Npix=60; %using the default for now, but keeping as variable in case we want to change later

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

if nargin < 6 || isempty(thresh)
    thresh = 0.3;
end

if nargin < 7 || isempty(num_lags)
    num_lags = 6;
end

if nargin < 8 || isempty(save_vars)
    save_vars.to_save=0;
end

tic
stim_shift=permute(single(data.stim),[4 1 2 3])./127;

if to_shift==1
    stim_shift = shift_stim( stim_shift, data.ETtrace, stim_deltas );
    shiftstr = 'with shifts';
elseif to_shift==0
    stim_shift = shift_stim( stim_shift, zeros(size(data.ETtrace)), stim_deltas );
    shiftstr = 'no shifts';
else
    shiftstr = 'no shifts';
end

stim_shift_L=squeeze(stim_shift(:,:,:,1)*[-0.040] + stim_shift(:,:,:,2)*[0.5220]); % from 01_2022 calib
stim_shift_M=squeeze(stim_shift(:,:,:,1)*[0.0351] + stim_shift(:,:,:,2)*[-0.4625]); % from 01_2022 calib

[fwdc_L, fwdc_M, fwdc_diff] = deal(zeros(max(targets),num_lags,Npix,Npix));

for pixX=1:Npix
    for pixY=1:Npix
        fwdc_inds_Lp{pixX,pixY} = find(stim_shift_L(use_inds,pixX,pixY)>thresh);
        fwdc_inds_Lm{pixX,pixY} = find(stim_shift_L(use_inds,pixX,pixY)<-thresh);
        fwdc_inds_Mp{pixX,pixY} = find(stim_shift_M(use_inds,pixX,pixY)>thresh);
        fwdc_inds_Mm{pixX,pixY} = find(stim_shift_M(use_inds,pixX,pixY)<-thresh);
    end
end
toc

disp('Now plotting!')
for cc=targets
	tic

	fwdcfig = figure;

	for curlag=1:num_lags
        for pixX=1:Npix
            for pixY=1:Npix
		        fwdc_L(cc,curlag,pixX,pixY) = mean(binned_SU(use_inds(fwdc_inds_Lp{pixX,pixY})+curlag,cc))-mean(binned_SU(use_inds(fwdc_inds_Lm{pixX,pixY})+curlag,cc));
		        fwdc_M(cc,curlag,pixX,pixY) = mean(binned_SU(use_inds(fwdc_inds_Mp{pixX,pixY})+curlag,cc))-mean(binned_SU(use_inds(fwdc_inds_Mm{pixX,pixY})+curlag,cc));
                fwdc_diff(cc,curlag,:,:) = squeeze(fwdc_L(cc,curlag,:,:))-squeeze(fwdc_M(cc,curlag,:,:));
            end
        end
	end
	curlim=max(abs(fwdc_L(:)')); 
        if curlim==0; curlim=0.1; end % avoids plotting bugs if a bad STA is included in a large set of plots

	for curlag=1:num_lags

		subplot(num_lags,3,1+(curlag-1)*3)
		imagesc(squeeze(squeeze(fwdc_L(cc,curlag,:,:)))); clim([-curlim curlim]); pbaspect([1 1 1])
        if curlag == 1; title('L'); end
        ylabel(['Lag ' num2str(curlag)]);

		subplot(num_lags,3,2+(curlag-1)*3)
		imagesc(squeeze(squeeze(fwdc_M(cc,curlag,:,:)))); clim([-curlim curlim]); pbaspect([1 1 1])
        if curlag == 1; title('M'); end

 		subplot(num_lags,3,3+(curlag-1)*3)
		imagesc(squeeze(squeeze(fwdc_diff(cc,curlag,:,:)))); clim([-curlim curlim]); pbaspect([1 1 1])
        if curlag == 1; title('L-M'); end

		if curlag == 1
			xlabel('Lum          L-M          S'); 
            if cc<=nSU;
                titlestr = sprintf('SU %0d: ch# %d; %s', cc, spk_ch(cc), shiftstr);
            else
                titlestr = sprintf('MU %0d: ch# %d; %s', cc-nSU, spk_ch(cc), shiftstr);
            end
            sgtitle(titlestr)
        end
		
		colormap(jet); 

	end
	toc
	%    pause
    if save_vars.to_save==1
        % save the STA as pdf saveas(figure, filename)
    	saveas(fwdcfig,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_STA.pdf'])
        % save the STA as a jpg
        % saveas(STA,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_STA.jpg'])
    end
end

