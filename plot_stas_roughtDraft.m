function plot_stas_roughtDraft(STA, Robs, chans, clusters, mainKofiko_fname, trialTypeOfInterest)

if nargin ==6
    [a,filenameP,~] = fileparts(mainKofiko_fname);
    kofiko_subfolder = fullfile(a,filenameP);
    trial = loadKofikoTrialData(kofiko_subfolder,mainKofiko_fname,filenameP);

    StimRect = cat(1,trial.m_aiStimulusRect);

    if strcmpi(trialTypeOfInterest, 'ccloud')

        isTrialOfInterest = strcmpi({trial.m_strTrialType}', 'Dual Stim') &  [trial.DualstimPrimaryuseRGBCloud]' == 8;
        shift = [30,30];
    elseif strcmpi(trialTypeOfInterest, 'hartley')
        isTrialOfInterest = strcmpi({trial.m_strTrialType}', 'Dual Stim') &  [trial.DualstimPrimaryuseRGBCloud]' == 6;
        shift = [0,0];
    end

    modalStimRect = mode(StimRect(isTrialOfInterest,:));
  
    %modalStimRect = mode(StimRect);


    x_left = modalStimRect(1) - 960;
    y_top = modalStimRect(2) - 540;
    w = modalStimRect(3) - modalStimRect(1);
    h = modalStimRect(4) - modalStimRect(2);

    x_right  = x_left + w;
    y_bottom = y_top + h;

else
    shift = [30,30];
    x_left = 1; x_right = 60;
    y_top =1; y_bottom = 60;
end

xt = linspace(x_left, x_right, 5);
yt = linspace(y_top,  y_bottom, 5);
xticklabs = cellstr(num2str(round(xt(:))));
yticklabs = cellstr(num2str(round(yt(:))));

chrom_chans = {'Lum.', 'L-M', 'S'};

lags = 2:7;
for i = 1:size(STA,1)

    figure;
    subplot_counter = 1;
    %tl = tiledlayout(3,6,'TileSpacing','tight','Padding','tight');
    chanNum =chans(i) - 1;
    clusterNum =clusters(i);

    N = sum(Robs(i,:));
    
    figTitle = sprintf('Chan %d, unit %d, spikes: %d', ...
        chanNum, clusterNum, N);
    
    set(gcf, 'Name', figTitle)
   
    for c = 1:3
        for j = lags
            %ax = nexttile;
             ax = subplot(3,6,subplot_counter);

            A = circshift(squeeze(STA(i,:,:,c,j+1)), shift);
            imagesc([x_left x_right], [y_top y_bottom], A, 'Parent', ax);
            axis(ax,'square')
            colormap(ax,gray)
            set(ax,'FontSize',16, ...
                'TickDir','out', ...
                'TickLength',[0.03 0.03], ...
                'XTick',xt, ...
                'YTick',yt, ...
                'LineWidth', 2)

            if c == 1
                title(ax,['Lag ' num2str(1*(j))])
            end

            if j == lags(1)
                ylabel(ax,chrom_chans{c})
            end

            if c == 1 && j == lags(1)
                set(ax,'XTickLabel',xticklabs, ...
                    'YTickLabel',yticklabs)
            else
                set(ax,'XTickLabel',[], ...
                    'YTickLabel',[])
            end
            subplot_counter = subplot_counter+1;
        end
    end



    %figname = ['Channel' num2str(chanNum) 'Unit' num2str(clusterNum)];
    %exportgraphics(gcf, fullfile('/Users/greenemj/Sprout/260728',[figname '.png']), 'Resolution',300)

end
end