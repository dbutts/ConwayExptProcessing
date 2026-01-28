function stas = get_sta(data, ExptInfo, targets, use_inds, to_shift, stim_deltas, num_lags, save_vars)
% Usage: stas = get_sta(data, ExptInfo, <targets>, <use_inds>, <to_shift>, <stim_deltas>, <num_lags>, <save_vars>)

% spikes & metadata
binned_SU = [single(data.Robs'), single(data.RobsMU')];
spk_ch    = [data.Robs_probe_ID; data.RobsMU_probe_ID];
nSU       = size(data.Robs,1);
nMU       = size(data.RobsMU,1);
nTOT      = nSU + nMU;

% Phy IDs from ExptInfo
unitIDs_phy = [ExptInfo.spk_ID_SU(:); ExptInfo.spk_ID_MU(:)];
isSU = [true(nSU,1); false(nMU,1)];

% targets (accept internal indices or Phy IDs)
if nargin < 3 || isempty(targets)
    targets = 1:nTOT;
else
    want_as_phy = any(targets > nTOT) || any(~ismember(targets, 1:nTOT));
    if want_as_phy && all(ismember(targets, unitIDs_phy))
        [~, targets] = ismember(targets(:), unitIDs_phy);
        targets = targets(:).';
    end
end

% indices
if nargin < 4 || isempty(use_inds)
    use_inds = data.valid_data;
    if numel(use_inds) > 10, use_inds(end-9:end) = []; end
end

% shifts / params (defaults)
if nargin < 5 || isempty(to_shift), to_shift = 0; end
if nargin < 6 || isempty(stim_deltas)
    NTet = size(data.ETtrace,2);
    stim_deltas = zeros(2, NTet);
end
if nargin < 7 || isempty(num_lags), num_lags = 6; end
if nargin < 8 || isempty(save_vars), save_vars.to_save = 0; end
if ~isfield(save_vars,'outputdir') || isempty(save_vars.outputdir), save_vars.outputdir = pwd; end
if ~isfield(save_vars,'titlestr'), save_vars.titlestr = ''; end
if ~isfield(save_vars,'save_mat'), save_vars.save_mat = 0; end
if ~isfield(save_vars,'close_figs'), save_vars.close_figs = false; end

% stimulus
stim_shift = permute(data.stim,[4 1 2 3]);
if to_shift==1
    stim_shift = shift_stim(stim_shift, data.ETtrace, stim_deltas);
    shiftstr = 'with shifts';
else
    stim_shift = shift_stim(stim_shift, zeros(size(data.ETtrace)), stim_deltas);
    shiftstr = 'no shifts';
end
stim2 = single(reshape(stim_shift, size(stim_shift,1), 3*60*60))./127;

% guard indices and precompute stimulus rows used
NT = size(binned_SU,1);
use_inds = use_inds(:);
use_inds = use_inds(use_inds >= 1 & use_inds <= NT - num_lags);
S = stim2(use_inds,:);

% SU/MU ordinals
idx_SU = find(isSU); idx_MU = find(~isSU);
ord_SU = zeros(nTOT,1); ord_SU(idx_SU) = 1:numel(idx_SU);
ord_MU = zeros(nTOT,1); ord_MU(idx_MU) = 1:numel(idx_MU);

% output
stas = zeros(nTOT, num_lags, 3*60*60, 'single');

exptname = 'rec';
if isfield(ExptInfo,'exptname') && ~isempty(ExptInfo.exptname)
    exptname = ExptInfo.exptname;
end

disp('Now plotting!')
for cc = targets
    cur_STA1 = zeros(num_lags, size(S,2), 'single');

    if isSU(cc)
        nameTag = sprintf('STA SU %03d (Phy %d)', ord_SU(cc), unitIDs_phy(cc));
    else
        nameTag = sprintf('STA MU %03d (Phy %d)', ord_MU(cc), unitIDs_phy(cc));
    end
    hFig = figure('Name', nameTag, 'Visible','on');

    for curlag = 1:num_lags
        cur_STA1(curlag,:) = single(binned_SU(use_inds+curlag, cc))' * S;

        curlim = max(abs(cur_STA1(1:curlag,:)), [], 'all');
        if curlim == 0, curlim = 0.1; end

        cur_StA2 = reshape(cur_STA1(curlag,:),60,180);
        cur_StA2(:,[60 120]) = curlim;

        subplot(6,3,1+(curlag-1)*3)
        imagesc(cur_StA2); caxis([-curlim curlim]); pbaspect([3 1 1])
        if curlag == 1, xlabel('Lum  L-M  S'); end

        cur_StA3  = reshape(cur_STA1(curlag,:),60*60,3);
        cur_StA_L = cur_StA3*[-0.040, 0.5220, 0]';
        cur_StA_M = cur_StA3*[ 0.0351,-0.4625, 0]';
        cur_StA_L2 = reshape(cur_StA_L,60,60);
        cur_StA_M2 = reshape(cur_StA_M,60,60);

        subplot(6,3,2+(curlag-1)*3)
        imagesc(cur_StA_L2); caxis([-curlim curlim]*0.522); pbaspect([1 1 1])
        if curlag == 1, xlabel('L'); end

        subplot(6,3,3+(curlag-1)*3)
        imagesc(cur_StA_M2); caxis([-curlim curlim]*0.4625); pbaspect([1 1 1])
        xlabel(['Lag ' num2str(curlag)])
        if curlag == 1, xlabel('M'); end

        if curlag == 1
            if isSU(cc)
                titlestr = sprintf('SU %d (Phy %d): ch %d; %s', ord_SU(cc), unitIDs_phy(cc), spk_ch(cc), shiftstr);
            else
                titlestr = sprintf('MU %d (Phy %d): ch %d; %s', ord_MU(cc), unitIDs_phy(cc), spk_ch(cc), shiftstr);
            end
            title(titlestr)
        end

        colormap(gray);
        drawnow limitrate
    end

    stas(cc,:,:) = cur_STA1;

    % filenames
    if isSU(cc), utype = 'SU'; uidx = ord_SU(cc); else, utype = 'MU'; uidx = ord_MU(cc); end
    base = sprintf('%s_%s_phy%05d_%s%03d_STA', save_vars.titlestr, exptname, unitIDs_phy(cc), utype, uidx);

    % ensure output dir if saving
    if (save_vars.to_save == 1 || save_vars.save_mat == 1) && ~exist(save_vars.outputdir,'dir')
        mkdir(save_vars.outputdir);
    end

    % save PDF
    if data.to_save == 1
        saveas(hFig, fullfile(data.outputdir, [exptname, '_PhyClus', num2str(unitIDs_phy(cc)), utype, num2str(uidx), '.pdf']));
    end

    % save per-unit MAT
    if save_vars.save_mat == 1
        matname = fullfile(save_vars.outputdir, [base '.mat']);
        STA = struct();
        STA.version        = 1;
        STA.exptname       = exptname;
        STA.titlestr       = save_vars.titlestr;
        STA.phy_id         = unitIDs_phy(cc);
        STA.unit_type      = utype;
        STA.unit_ordinal   = uidx;
        STA.array_index    = cc;
        STA.channel        = spk_ch(cc);
        STA.isSU           = isSU(cc);
        STA.num_lags       = num_lags;
        STA.grid_sz        = [60 60];
        STA.color_channels = {'R','G','B'};
        STA.transform_L    = [-0.040,  0.5220,  0];
        STA.transform_M    = [ 0.0351, -0.4625, 0];
        STA.to_shift       = logical(to_shift);
        STA.shift_desc     = shiftstr;
        STA.stim_deltas    = stim_deltas;
        STA.use_inds       = use_inds;
        STA.NT             = NT;
        STA.sta_flat       = cur_STA1;        % [num_lags x 10800] single
        STA.display_name   = nameTag;
        STA.saved_at       = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        try
            save(matname, 'STA');
        catch
            save(matname, 'STA', '-v7.3');
        end
    end

    if save_vars.close_figs
        close(hFig);
    end
end