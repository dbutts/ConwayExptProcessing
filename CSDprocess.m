function [csd, lfp] = CSDprocess( ExptRecording, trial_select, save_images, exptname, savedir, verbose )
%
% Usage: [csd, lfp] = CSDprocess( ExptRecording, <trial_select>, <save_images>, <exptname>, <savedir>, <verbose> )
%
% Trial select could be range or single number. If single number, positive
% will be 1:num and negative will be the last num trials

if nargin < 6
	verbose = 1;
end
if (nargin < 5) || isempty(savedir)
	savedir = '';
end
if (nargin < 4) || isempty(exptname)
	exptname = '';
	if save_images == 1
		disp( 'No experiment name specified: not saving images.')
		save_images = 0;
	end
else
	% Scrub expt-name
	for ii = 1:length(exptname)
		if exptname(ii) == '_'
			exptname(ii) = '-';
		end
	end
end
if nargin < 3
	save_images=0;
end

nChans=24;

%% Now we process LFPs
targ_trials=[];
for tt=1:size(ExptRecording,1)
	if strcmp(ExptRecording{tt, 1}.m_strTrialType, 'Dual Stim') && ...
			ExptRecording{tt, 1}.m_bMonkeyFixated==1
		targ_trials=[targ_trials,tt];    
	end
end
Ntr = length(targ_trials);
if verbose
	fprintf( '%d valid trials identified\n', Ntr)
end

if nargin > 1
	% then there is a trial select if length > 0
	if length(trial_select) == 1
		if trial_select > 0
			if trial_select >= Ntr
				if verbose
					disp( 'Using all trials.' )
				end
			else
				targ_trials = targ_trials(1:trial_select);
			end
		else
			if abs(trial_select) >= Ntr
				if verbose
					disp( 'Using all trials.' )
				end
			else
				targ_trials = targ_trials(Ntr+trial_select:end);
			end
		end
	elseif length(trial_select) == 2
		targ_trials = targ_trials(trial_select(1):trial_select(2));
	elseif length(trial_select) > 2
		targ_trials = targ_trials(trial_select);
	end
end

%% Get CSD
all_lfps=zeros(nChans,800,length(targ_trials));

iii=1;

for tt=targ_trials
	%if ExperimentRecording{tt, 1}.CSDtrigframe==1
	if ~isempty(ExptRecording{tt,4})   % BANDAID to get rid of bad first trial in ExperimentRecord
		all_lfps(:,:,iii)=ExptRecording{tt,4}(:,1:800);
		all_lfptimes(iii)=ExptRecording{tt,2};
		iii=iii+1;        
	end
end


%%
cur_lfp_inds=[1:iii-1];
%cur_lfp_inds=[160:iii-1];
%twin = all_lfptimes([cur_lfp_inds(1) cur_lfp_inds(end)])/60;

%%
vars.Fs = 1000;
vars.BrainBound = 1;
vars.ChanSep = 0.05;
vars.diam = 2; %0.5
CSD = PettersenCSD(all_lfps(:,:,cur_lfp_inds),'spline',vars);

%%
expt_CSDs = mean(CSD,3);
CSDfig=figure;
imagesc(expt_CSDs)
title([exptname ' iCSD'])
ylabel('Probe number'); xlabel('Time (ms)')
xlim([0 550]); colorbar
ax=gca; ax.XTickLabel={'-100','0','100','200','300','400','500'};
%ax=gca; ax.XTickLabel={'0','100','200','300','400','500','600'};
%caxis([-4e5 1e5])
%ax=gca; ax.XTickLabel={'-50','50','150','250','350','450'};
%colormap(hotcold); caxis([-max(expt_CSDs(:))/3 max(expt_CSDs(:))/3]);

if save_images
	saveas(CSDfig,[savedir exptname '_iCSD.pdf'])
	saveas(CSDfig,[savedir exptname '_iCSD.png'])
end

CSDfig2=figure;
imagesc(expt_CSDs)
title([exptname ' iCSD'])
ylabel('Probe number'); xlabel('Time (ms)')
xlim([0 250]); colorbar
ax=gca; ax.XTickLabel={'0','50','100','150','200','250'};

if save_images
	saveas(CSDfig2,[savedir exptname '_iCSD_zoom.pdf'])
	saveas(CSDfig2,[savedir exptname '_iCSD_zoom.png'])
end

%%
lfps2=squeeze(nanmean(all_lfps(:,:,cur_lfp_inds), 3));
LFPfig=figure;
imagesc(lfps2)
title([exptname ' LFP'])
ylabel('Probe number'); xlabel('Time (ms)')
xlim([0 550]); colorbar
ax=gca; ax.XTickLabel={'-100','0','100','200','300','400','500'};
%colormap(hotcold); caxis([-max(lfps2(:))/3 max(lfps2(:))/3]);

if save_images
	saveas(LFPfig,[savedir exptname '_LFP.pdf'])
	saveas(LFPfig,[savedir exptname '_LFP.png'])
end
