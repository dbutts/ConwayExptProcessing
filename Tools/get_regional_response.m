function get_regional_response(data, stas, target_SUs, apply_ETshifts, num_lags, save_vars, type1, parameters1, type2, parameters2, type3, parameters3, type4, parameters4)
%
% Usage: stas = get_sta( data, <targets>, <apply_ETshifts>, <stim_deltas> )
%
% PLAN COMING BACK:
% -- Make STA function that uses the shifted stim and spikes (and valid data)
%    It then calls the STA function, as well as shifting function, and reshaping and so on 

%% Now we process the data to make STAs
binned_SU = [single(data.Robs'), single(data.RobsMU')];
spk_ch = [data.Robs_probe_ID; data.RobsMU_probe_ID];
nSU=size(data.Robs,1); 

% 
% if nargin < 2 
%     target_SUs = 1:size(binned_SU,2); % default to plotting all SUs and MUs
% end
% 
% if nargin <3 || isempty(use_inds)
%     use_inds = data.valid_data;
%     use_inds(end-10:end) = []; %cut last few indices to avoid artifacts
% end
% 
% if nargin < 4 || isempty(apply_ETshifts)
%     apply_ETshifts = 0; 
% end
% 
% if nargin < 5
%     NT=size(data.ETtrace,2);
%     stim_deltas = zeros(2,NT); 
% end
% 
% if nargin < 6 || isempty(num_lags)
%     num_lags = 6;
% end
% 
% if nargin < 7 || isempty(save_vars)
%     save_vars.to_save=0;
% end
% 
% stim_shift=permute(data.stim,[4 1 2 3]);
% 
if apply_ETshifts==1
    % stim_shift = shift_stim( stim_shift, data.ETtrace, stim_deltas );
    shiftstr = 'with shifts';
elseif apply_ETshifts==0
    % stim_shift = shift_stim( stim_shift, zeros(size(data.ETtrace)), stim_deltas );
    shiftstr = 'no shifts';
else
    shiftstr = 'no shifts';
end
% 
% %reshaping stimuli to get STAs with matrix multiplication, which is faster
% stim2 = single(reshape(stim_shift,size(stim_shift,1),3*60*60))./127;
% 
% stas=zeros(max(target_SUs),num_lags,3*60*60);

area1 = zeros(5,4);
area2 = zeros(5,4);
area3 = zeros(5,4);
area4 = zeros(5,4);
% rec1_LMy = [];
% rec1_Sy = [];
% rec1_My = [];
% rec1_Ly = [];
% 
% rec2_lumy = [];
% rec2_LMy = [];
% rec2_Sy = [];
% rec2_My = [];
% rec2_Ly = [];
% 
% rec3_lumy = [];
% rec3_LMy = [];
% rec3_Sy = [];
% rec3_My = [];
% rec3_Ly = [];
% 
% rec4_lumy = [];
% rec4_LMy = [];
% rec4_Sy = [];
% rec4_My = [];
% rec4_Ly = [];

disp('Now plotting!')
for cc=target_SUs
	tic
    % % cameron added 2 lines: 
    % indexed_spikes = binned_SU(use_inds,cc);
    % total_spikes = sum(indexed_spikes == 1);
    
	STA = figure;
	% for curlag=1:num_lags
	% 	cur_STA1(curlag,:) = binned_SU(use_inds+curlag,cc)' * stim2(use_inds,:);
	% end

    % stas(cc,:, :)=cur_STA1;
    cur_STA1 = stas(cc,:,:);
	curlim=max(abs(cur_STA1(:)')); 
        if curlim==0; curlim=0.1; end % avoids plotting bugs if a bad STA is included in a large set of plots

	for curlag=1:num_lags
        cur_STA1 = squeeze(cur_STA1);
		cur_STA2=reshape(cur_STA1(curlag,:),60,180);
		%curlim=150;%
		cur_STA2(:,[60 120])=curlim; % plots vertical lines dividing up the plots

		subplot(6,3,1+(curlag-1)*3)
		imagesc(cur_STA2); clim([-curlim curlim]); pbaspect([3 1 1])

        % hold on 
        % 
        % lumy1 = [lumy1, avg_rectangle(rectangle1, cur_STA2, 'r')];
        % rec2_lumy = [rec2_lumy, avg_rectangle(rectangle2, cur_STA2, 'g')];
        % rec3_lumy = [rec3_lumy, avg_rectangle(rectangle3, cur_STA2, 'b')];
        % rec4_lumy = [rec4_lumy, avg_rectangle(rectangle4, cur_STA2, 'c')];
        % 
        % rectangle1(1) = rectangle1(1) +60;
        % rectangle2(1) = rectangle2(1) +60;
        % rectangle3(1) = rectangle3(1) +60;
        % rectangle4(1) = rectangle4(1) +60;
        % 
        % rec1_LMy = [rec1_LMy, avg_rectangle(rectangle1, cur_STA2, 'r')];
        % rec2_LMy = [rec2_LMy, avg_rectangle(rectangle2, cur_STA2, 'g')];
        % rec3_LMy = [rec3_LMy, avg_rectangle(rectangle3, cur_STA2, 'b')];
        % rec4_LMy = [rec4_LMy, avg_rectangle(rectangle4, cur_STA2, 'c')];
        % 
        % rectangle1(1) = rectangle1(1) +60;
        % rectangle2(1) = rectangle2(1) +60;
        % rectangle3(1) = rectangle3(1) +60;
        % rectangle4(1) = rectangle4(1) +60;
        % 
        % rec1_Sy = [rec1_Sy, avg_rectangle(rectangle1, cur_STA2, 'r')];
        % rec2_Sy = [rec2_Sy, avg_rectangle(rectangle2, cur_STA2, 'g')];
        % rec3_Sy = [rec3_Sy, avg_rectangle(rectangle3, cur_STA2, 'b')];
        % rec4_Sy = [rec4_Sy, avg_rectangle(rectangle4, cur_STA2, 'c')];
        % 
        % rectangle1(1) = rectangle1(1) -120;
        % rectangle2(1) = rectangle2(1) -120;
        % rectangle3(1) = rectangle3(1) -120;
        % rectangle4(1) = rectangle4(1) -120;
        
        % lumy = avg_circle(cur_STA2, circle1);
        area1(1,curlag) = avg_region(cur_STA2, type1, parameters1, 0); %lum
        area1(2,curlag) = avg_region(cur_STA2, type1, parameters1, 60); %LM
        area1(3,curlag) = avg_region(cur_STA2, type1, parameters1, 120); %S

        area2(1,curlag) = avg_region(cur_STA2, type2, parameters2, 0); %lum
        area2(2,curlag) = avg_region(cur_STA2, type2, parameters2, 60); %LM
        area2(3,curlag) = avg_region(cur_STA2, type2, parameters2, 120); %S

        area3(1,curlag) = avg_region(cur_STA2, type3, parameters3, 0); %lum
        area3(2,curlag) = avg_region(cur_STA2, type3, parameters3, 60); %LM
        area3(3,curlag) = avg_region(cur_STA2, type3, parameters3, 120); %S

        area4(1,curlag) = avg_region(cur_STA2, type4, parameters4, 0); %lum
        area4(2,curlag) = avg_region(cur_STA2, type4, parameters4, 60); %LM
        area4(3,curlag) = avg_region(cur_STA2, type4, parameters4, 120); %S

        cur_StA3=reshape(cur_STA1(curlag,:),60*60,3);
        cur_StA_L = cur_StA3*[-0.040, 0.5220,0]'; % from 01_2022 calib
		cur_STA_L2=reshape(cur_StA_L,60,60);
        cur_StA_M = cur_StA3*[0.0351, -0.4625, 0]'; % from 01_2022 calib
		cur_STA_M2=reshape(cur_StA_M,60,60);

        if curlag == 1
			xlabel('Lum          L-M          S'); 
            if cc<=nSU
                titlestr = sprintf('SU %0d: ch# %d; %s', cc, spk_ch(cc), shiftstr);
            else
                titlestr = sprintf('MU %0d: ch# %d; %s', cc-nSU, spk_ch(cc), shiftstr);
            end
			%title(['SU ' num2str(cc) ': ch#', num2str(spk_ch(cc)), 'spkID' num2str(spk_ID(cc))]); 
            title(titlestr)
            % % cameron added 1 line:
            % subtitle("Total exp time: " + length(use_inds)/240 + " seconds" + newline + "end" ...
            %     + newline + "total spikes: " + total_spikes);
        end

        subplot(6,3,2+(curlag-1)*3)
		imagesc(cur_STA_L2); clim([-curlim curlim].*.522); pbaspect([1 1 1]); % color limit based on L-isolating max

        % rec1_Ly = [rec1_Ly, avg_rectangle(rectangle1, cur_STA_L2, 'r')];
        % rec2_Ly = [rec2_Ly, avg_rectangle(rectangle2, cur_STA_L2, 'g')];
        % rec3_Ly = [rec3_Ly, avg_rectangle(rectangle3, cur_STA_L2, 'b')];
        % rec4_Ly = [rec4_Ly, avg_rectangle(rectangle4, cur_STA_L2, 'c')];
        area1(4,curlag) = avg_region(cur_STA_L2, type1, parameters1, 0); %L
        area2(4,curlag) = avg_region(cur_STA_L2, type2, parameters2, 0); %L
        area3(4,curlag) = avg_region(cur_STA_L2, type3, parameters3, 0); %L
        area4(4,curlag) = avg_region(cur_STA_L2, type4, parameters4, 0); %L

		if curlag == 1; xlabel('L'); end

		subplot(6,3,3+(curlag-1)*3)
		imagesc(cur_STA_M2); clim([-curlim curlim].*.4625); pbaspect([1 1 1])

        % rec1_My = [rec1_My, avg_rectangle(rectangle1, cur_STA_M2, 'r')];
        % rec2_My = [rec2_My, avg_rectangle(rectangle2, cur_STA_M2, 'g')];
        % rec3_My = [rec3_My, avg_rectangle(rectangle3, cur_STA_M2, 'b')];
        % rec4_My = [rec4_My, avg_rectangle(rectangle4, cur_STA_M2, 'c')];
        area1(5,curlag) = avg_region(cur_STA_M2, type1, parameters1, 0); %M
        area2(5,curlag) = avg_region(cur_STA_M2, type2, parameters2, 0); %M
        area3(5,curlag) = avg_region(cur_STA_M2, type3, parameters3, 0); %M
        area4(4,curlag) = avg_region(cur_STA_M2, type4, parameters4, 0); %M
		
		colormap(gray); xlabel(['Lag ' num2str(curlag)]);

        if curlag == 1; xlabel('M'); end


	end
	%    figtitle(['Probe ' num2str(Robs_probe_ID(cc)) '   Unit ' num2str(cc) ])
	toc
	%    pause
    if save_vars.to_save==1
        % save the STA as pdf saveas(figure, filename)
    	saveas(STA,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_plot.pdf'])
        % save the STA as a jpg
        % saveas(STA,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_STA.jpg'])
    end
end

activations = figure;
ylimits = [-300,500];

labels = ["Lum red area","L-M","S","L","M"];
% make_plots([lumy1, rec1_LMy, rec1_Sy, rec1_Ly, rec1_My], ydim, 1, labels)
make_plots(area1, area2, ylimits, 1, labels)

labels = ["Lum blue area", "L-M", "S", "L", "M"];
% make_plots([rec2_lumy, rec2_LMy, rec2_Sy, rec2_Ly, rec2_My], ydim, 6, labels)
make_plots(area2, area2, ylimits, 6, labels)

labels = ["Lum green area", "L-M", "S", "L", "M"];
% make_plots([rec3_lumy, rec3_LMy, rec3_Sy, rec3_Ly, rec3_My], ydim, 11, labels)
make_plots(area3, area2, ylimits, 11, labels)

labels = ["Lum yellow area", "L-M", "S", "L", "M"];
% make_plots([rec4_lumy, rec4_LMy, rec4_Sy, rec4_Ly, rec4_My], ydim, 16, labels)
make_plots(area4, area2, ylimits, 16, labels)

% plot(1:size(area1(1,:),2), area1(1,:))
% title('circle lum')

    if save_vars.to_save==1
        % save the STA as pdf saveas(figure, filename)
    	saveas(activations,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_activations.pdf'])
        % save the STA as a jpg
        % saveas(STA,[save_vars.outputdir save_vars.titlestr ' SU' num2str(cc) '_STA.jpg'])
    end


    function avg = avg_region(STA, type, parameters, xshift)
        if type(1) == "circle"
            avg = avg_circle(STA,type,parameters,xshift);
        elseif type(1) == "rectangle"
            avg = avg_rectangle(STA,type,parameters,xshift);
        elseif type(1) == "polygon"
            avg = avg_polygon(STA,type,parameters,xshift);
        end
    end

    function avg = avg_polygon(STA,type,parameters, xshift)
        xv = parameters(1,:) + xshift;
        yv = parameters(2,:);
        [X,Y] = meshgrid(1:size(STA,2),1:size(STA,1));
        xq = X(:);
        yq = Y(:);
        [in,~] = inpolygon(xq,yq,xv,yv);
        color = type(2);
        plot(xq(in), yq(in),'.','Color',color,'Markersize',1)    
        avg = mean(STA(in));
    end

    function avg = avg_circle(STA,type,parameters,xshift)
        xc = parameters(1) + xshift;
        yc = parameters(2);
        r = parameters(3);
        xv = linspace(1, size(STA,1), size(STA, 1));
        yv = linspace(1, size(STA,2), size(STA, 2));
        [X,Y] = ndgrid(xv, yv);
        Isoutside = find((X(:)-yc).^2+(Y(:)-xc).^2>r^2);
        STA(Isoutside) = nan;
        avg = nanmean(STA, "all");
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + xc;
        yunit = r * sin(th) + yc;
        hold on
        color = type(2);
        plot(xunit, yunit, color, 'LineWidth',1);
    end

    function avg = avg_rectangle(STA,type,parameters,xshift)
        color = type(2);
        parameters(1) = parameters(1) + xshift;
        avg = mean(STA(parameters(2):parameters(2)+parameters(4)-1, parameters(1):parameters(1)+parameters(3)-1), 'all');
        rectangle('Position', parameters, 'EdgeColor', color)
    end

    % function make_plots(list, ylimits, startingplot, labels)
    %     for i = 1:5
    %         subplot(4,5,startingplot+i-1)
    %         y = list((i-1)*6+1:(i-1)*6+6);
    %         plot(1:size(y, 2), y)
    %         ylim(ylimits);
    %         title(labels(i));
    %         ylabel('STA value');
    %         xlabel('lag number'); 
    %     end
    %     sgtitle('Activation in each rectangle')
    % end

    function make_plots(area, area2, ylimits, startingplot, labels)
        for i = 1:5
            subplot(4,5,startingplot+i-1)
            plot(1:size(area, 2), area(i,:)-area2(i,:))
            ylim(ylimits);
            title(labels(i));
            ylabel('STA value');
            xlabel('lag number'); 
            xlim([1,5])
        end
        sgtitle('Blue normalized activation in each region')
    end

end
 