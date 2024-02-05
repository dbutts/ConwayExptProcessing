function shifted_stim = shift_stim( stim, ETtrace, stim_deltas )
%
% Usage: shifted_stim = shift_stim( stim, ETtrace, <stim_deltas> )
% 
% stim passed in should be size [NT NX NY nlags]
% Note that shifts larger than 30 in any direction will be misinterpreted (circ)

NT = size(stim, 1);
shifted_stim = stim;

if nargin < 3
	stim_deltas = zeros(NT,2);
end

shifts = mod( round(-ETtrace)+30 + stim_deltas, 60);
shifts(isnan(shifts)) = 30;

tic
for i = 1:NT
	shifted_stim(i,:,:,:) = circshift( shifted_stim(i,:,:,:), shifts(1,i), 2 );
	shifted_stim(i,:,:,:) = circshift( shifted_stim(i,:,:,:), shifts(2,i), 3 );
end
toc

% stim2 = single(reshape(stim_shift,size(stim_shift,1),3*60*60))./127;
