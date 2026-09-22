function plot_stas_roughtDraft(STA, Robs,chans, clusters, mainKofiko_fname)

if nargin ==5
    [a,filenameP,~] = fileparts(mainKofiko_fname);
    kofiko_subfolder = fullfile(a,filenameP);
    trial = loadKofikoTrialData(kofiko_subfolder,mainKofiko_fname,filenameP);
    modalStimRect = mode(cat(1,trial.m_aiStimulusRect));

    x_left = modalStimRect(1) - 960;
    y_top = modalStimRect(2) - 540;
    w = modalStimRect(3) - modalStimRect(1);
    h = modalStimRect(4) - modalStimRect(2);

    x_right  = x_left + w;
    y_bottom = y_top + h;

else
    x_left = 1; x_right = 60;
    y_top =1; y_bottom = 60;
end

xt = linspace(x_left, x_right, 5);
yt = linspace(y_top,  y_bottom, 5);
xticklabs = cellstr(num2str(round(xt(:))));
yticklabs = cellstr(num2str(round(yt(:))));

chrom_chans = {'Lum.', 'L-M', 'S'};
for i = 1:size(STA.DKL,1)

    figure
    tl = tiledlayout(3,6,'TileSpacing','tight','Padding','tight');

    for c = 1:3
        for j = 3:8
            ax = nexttile;

            A = circshift(squeeze(STA.DKL(i,:,:,c,j)), [30 30]);
            imagesc(ax, [x_left x_right], [y_top y_bottom], A);
            axis(ax,'square')
            colormap(ax,gray)
            set(ax,'FontSize',16, ...
                'TickDir','out', ...
                'TickLength',[0.03 0.03], ...
                'XTick',xt, ...
                'YTick',yt, ...
                'LineWidth', 2)

            if c == 1
                title(ax,['Lag ' num2str(1*(j-1))])
            end

            if j == 3
                ylabel(ax,chrom_chans{c})
            end

            if c == 1 && j == 3
                set(ax,'XTickLabel',xticklabs, ...
                    'YTickLabel',yticklabs)
            else
                set(ax,'XTickLabel',[], ...
                    'YTickLabel',[])
            end
        end
    end

    chanNum =chans(i) -1;
    clusterNum =clusters(i);

    N = sum(Robs(i,:));
    sgtitle(tl, sprintf('Chan %d, unit %d, spikes: %d', ...
        chanNum, clusterNum, N));

    %figname = ['Channel' num2str(chanNum) 'Unit' num2str(clusterNum)];
    %exportgraphics(gcf, fullfile('/Users/greenemj/Sprout/260728',[figname '.png']), 'Resolution',300)

end
end