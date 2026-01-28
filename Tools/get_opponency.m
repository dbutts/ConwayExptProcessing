function stas = get_opponency( data, targets, use_inds, to_shift, stim_deltas, num_lags, save_vars, rectangle1, rectangle2, rectangle3, rectangle4)
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

rec1_lumy = [];
rec1_LMy = [];
rec1_Sy = [];
rec1_My = [];
rec1_Ly = [];

rec2_lumy = [];
rec2_LMy = [];
rec2_Sy = [];
rec2_My = [];
rec2_Ly = [];

rec3_lumy = [];
rec3_LMy = [];
rec3_Sy = [];
rec3_My = [];
rec3_Ly = [];

rec4_lumy = [];
rec4_LMy = [];
rec4_Sy = [];
rec4_My = [];
rec4_Ly = [];

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
		cur_STA2=reshape(cur_STA1(curlag,:),60,180);
		%curlim=150;%
		cur_STA2(:,[60 120])=curlim; % plots vertical lines dividing up the plots

		subplot(6,3,1+(curlag-1)*3)
		imagesc(cur_STA2); clim([-curlim curlim]); pbaspect([3 1 1])

        hold on 
        
        rec1_lumy = [rec1_lumy, make_rec_list(rectangle1, cur_STA2, 'r')];
        rec2_lumy = [rec2_lumy, make_rec_list(rectangle2, cur_STA2, 'g')];
        rec3_lumy = [rec3_lumy, make_rec_list(rectangle3, cur_STA2, 'b')];
        rec4_lumy = [rec4_lumy, make_rec_list(rectangle4, cur_STA2, 'c')];

        rectangle1(1) = rectangle1(1) +60;
        rectangle2(1) = rectangle2(1) +60;
        rectangle3(1) = rectangle3(1) +60;
        rectangle4(1) = rectangle4(1) +60;

        rec1_LMy = [rec1_LMy, make_rec_list(rectangle1, cur_STA2, 'r')];
        rec2_LMy = [rec2_LMy, make_rec_list(rectangle2, cur_STA2, 'g')];
        rec3_LMy = [rec3_LMy, make_rec_list(rectangle3, cur_STA2, 'b')];
        rec4_LMy = [rec4_LMy, make_rec_list(rectangle4, cur_STA2, 'c')];

        rectangle1(1) = rectangle1(1) +60;
        rectangle2(1) = rectangle2(1) +60;
        rectangle3(1) = rectangle3(1) +60;
        rectangle4(1) = rectangle4(1) +60;

        rec1_Sy = [rec1_Sy, make_rec_list(rectangle1, cur_STA2, 'r')];
        rec2_Sy = [rec2_Sy, make_rec_list(rectangle2, cur_STA2, 'g')];
        rec3_Sy = [rec3_Sy, make_rec_list(rectangle3, cur_STA2, 'b')];
        rec4_Sy = [rec4_Sy, make_rec_list(rectangle4, cur_STA2, 'c')];

        rectangle1(1) = rectangle1(1) -120;
        rectangle2(1) = rectangle2(1) -120;
        rectangle3(1) = rectangle3(1) -120;
        rectangle4(1) = rectangle4(1) -120;

        cur_StA3=reshape(cur_STA1(curlag,:),60*60,3);
        cur_StA_L = cur_StA3*[-0.040, 0.5220,0]'; % from 01_2022 calib
		cur_STA_L2=reshape(cur_StA_L,60,60);
        cur_StA_M = cur_StA3*[0.0351, -0.4625, 0]'; % from 01_2022 calib
		cur_STA_M2=reshape(cur_StA_M,60,60);

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
		imagesc(cur_STA_L2); clim([-curlim curlim].*.522); pbaspect([1 1 1]); % color limit based on L-isolating max

        rec1_Ly = [rec1_Ly, make_rec_list(rectangle1, cur_STA_L2, 'r')];
        rec2_Ly = [rec2_Ly, make_rec_list(rectangle2, cur_STA_L2, 'g')];
        rec3_Ly = [rec3_Ly, make_rec_list(rectangle3, cur_STA_L2, 'b')];
        rec4_Ly = [rec4_Ly, make_rec_list(rectangle4, cur_STA_L2, 'c')];

		if curlag == 1; xlabel('L'); end

		subplot(6,3,3+(curlag-1)*3)
		imagesc(cur_STA_M2); clim([-curlim curlim].*.4625); pbaspect([1 1 1])

        rec1_My = [rec1_My, make_rec_list(rectangle1, cur_STA_M2, 'r')];
        rec2_My = [rec2_My, make_rec_list(rectangle2, cur_STA_M2, 'g')];
        rec3_My = [rec3_My, make_rec_list(rectangle3, cur_STA_M2, 'b')];
        rec4_My = [rec4_My, make_rec_list(rectangle4, cur_STA_M2, 'c')];

		
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

figure;

ydim = [-2000,1000];
labels = ["Lum - rectangle 1 (red)","L-M","S","L","M"];
make_plots([rec1_lumy, rec1_LMy, rec1_Sy, rec1_Ly, rec1_My], ydim, 1, labels)

labels = ["Lum - rectangle 2 (green)", "L-M", "S", "L", "M"];
make_plots([rec2_lumy, rec2_LMy, rec2_Sy, rec2_Ly, rec2_My], ydim, 6, labels)

labels = ["Lum - rectangle 3 (blue)", "L-M", "S", "L", "M"];
make_plots([rec3_lumy, rec3_LMy, rec3_Sy, rec3_Ly, rec3_My], ydim, 11, labels)

labels = ["Lum - rectangle 4 (cyan)", "L-M", "S", "L", "M"];
make_plots([rec4_lumy, rec4_LMy, rec4_Sy, rec4_Ly, rec4_My], ydim, 16, labels)


    function list = make_rec_list(rec_in, cur_STA2, color)
        list = mean(cur_STA2(rec_in(2):rec_in(2)+rec_in(4)-1, rec_in(1):rec_in(1)+rec_in(3)-1), 'all');
        rectangle('Position', rec_in, 'EdgeColor', color)
    end

    function make_plots(list, ylimits, startingplot, labels)
        for i = 1:5
            subplot(4,5,startingplot+i-1)
            y = list((i-1)*6+1:(i-1)*6+6);
            plot(1:size(y, 2), y)
            ylim(ylimits);
            title(labels(i));
            ylabel('STA value');
            xlabel('lag number'); 
        end
        sgtitle('Activation in each rectangle')
    end

end
 