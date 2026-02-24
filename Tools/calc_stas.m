function stas = calc_stas(data, cell_list, use_inds, to_shift, stim_deltas, num_lags, save_vars)
%function stas = get_sta2(data, ExptInfo, cell_list, use_inds, to_shift, stim_deltas, num_lags, save_vars)
% Usage: stas = get_sta2(data, ExptInfo, <cell_list>, <use_inds>, <to_shift>, <stim_deltas>, <num_lags>, <save_vars>)

% spikes & metadata
binned_SU = [single(data.Robs'), single(data.RobsMU')];
spk_ch    = [data.Robs_probe_ID; data.RobsMU_probe_ID];
nSU       = size(data.Robs,1);
nMU       = size(data.RobsMU,1);
nTOT      = nSU + nMU;
NC = length(cell_list);

% Phy IDs from ExptInfo
%unitIDs_phy = [ExptInfo.spk_ID_SU(:); ExptInfo.spk_ID_MU(:)];
isSU = [true(nSU,1); false(nMU,1)];

% cell_list (accept internal indices or Phy IDs)
if nargin < 2 || isempty(cell_list)
    cell_list = 1:nTOT;
else
    %want_as_phy = any(cell_list > nTOT) || any(~ismember(cell_list, 1:nTOT));
    %if want_as_phy && all(ismember(cell_list, unitIDs_phy))
    %    [~, cell_list] = ismember(cell_list(:), unitIDs_phy);
    cell_list = cell_list(:).';
    %end
end

% indices
if nargin < 3 || isempty(use_inds)
    use_inds = data.valid_data;
    if numel(use_inds) > 10, use_inds(end-9:end) = []; end
end

% shifts / params (defaults)
if nargin < 4 || isempty(to_shift), to_shift = 0; end
if nargin < 5 || isempty(stim_deltas)
    NTet = size(data.ETtrace,2);
    stim_deltas = zeros(2, NTet);
end
if nargin < 6 || isempty(num_lags), num_lags = 8; end
if nargin < 7 || isempty(save_vars), save_vars.to_save = 0; end
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
%idx_SU = find(isSU); idx_MU = find(~isSU);
%ord_SU = zeros(nTOT,1); ord_SU(idx_SU) = 1:numel(idx_SU);
%ord_MU = zeros(nTOT,1); ord_MU(idx_MU) = 1:numel(idx_MU);

% output
stas = zeros(NC, num_lags, 3*60*60, 'single');

%exptname = 'rec';
%if isfield(ExptInfo,'exptname') && ~isempty(ExptInfo.exptname)
%    exptname = ExptInfo.exptname;
%end

for ii = 1:length(cell_list)
		cc = cell_list(ii);
		fprintf('Cell %d\n', cc)
    for curlag = 1:num_lags
				%disp(curlag)
        stas(ii, curlag, :) = single(binned_SU(use_inds+curlag, cc))' * S;
		end
		stas(ii,:,:) = stas(ii,:,:) / max([1, sum(binned_SU(use_inds))]);
end

stas = reshape(stas, NC, num_lags, 60, 60, 3 ); 