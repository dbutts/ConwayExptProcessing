function [outmat, seedmat]=mk_ColorBarStim(pvec, width, stim_frames, nbars, seedmat_in)
% Example usage: 
%[outmat, seedmat]=GenerateTernaryStim(.35, 4, 0, [2 4], 10); imagesc(outmat(:,:)); colormap(gray(255))

% p: probability of black or white bars, should be below .5
% NOTE: All of these inputs should be integers
% width: width of the bars
% tilt: degree of orientation of the stimulus, in degrees clockwise
% stim_duration: number of frames the stimulus is repeated for
    % a single input means the number of frames is constant
    % range input (ex. [1 4]) randomly generates stimulus durations within that range
% run_duration: time in seconds you want the stimulus to run
    %If stimduration is >1 frame and randomly varied, run duration becaomes an expected value

% Outputs:
% outmat: height by bars by frames matrix
% seedmat: structure array containing metadata about the stimulus.
    % replicates the exact stimulus if input into this function

gray=127;     % value of gray on scale
%res=80;         % resolution in number of pixels across. res/width=number of bars
%midres=res/2;   % midpoint of resolution
%scaf=48;        % scaling factor for window; alternatively, use %scaf=round(res/7);
%ap1=midres-scaf/2; ap2=midres+scaf/2-1;
stim_duration= 1; %frames to be repeated - just want 1 unique frame for each refresh

if sum(pvec)~=1; pvec=pvec./sum(pvec); end
if nargin < 3; stim_frames=1; end
if nargin < 4; nbars=30; end

if nargin == 6 % Loads in information from seed matrix
pvec = seedmat_in.pvec;
width = seedmat_in.width;
%tilt = seedmat_in.tilt;
%stim_size = seedmat_in.stim_size;
end

iterations=stim_frames; %ceil(run_duration*60/mean(stim_duration)); %turns run duration from seconds into frames

randseed=rand(stim_frames, nbars);
barmat = zeros(stim_frames, nbars);
pvec_edges=[0 cumsum(pvec)];
for pp=1:7
    barmat(randseed>pvec_edges(pp) & randseed<pvec_edges(pp+1))=pp;
end
%outmat=repmat(barmat,1,1,nbars); expands stimuli to make a square texture
outmat=repmat(barmat,1,1,2);

seedmat.stims = squeeze(outmat(:,:,1));
seedmat.pvec = pvec;
seedmat.width = width;

% for i=1:iterations %looping through iterations; stimulus duration is added
%     if nargin == 5
%         durseed = seedmat_in.stimdurs(i,1);
%         barmat = repmat(seedmat_in.stims(i,:), res, 1);
%     else
%         randseed=rand(1,res);
% 
%         barmat=zeros(res, res); %preallocate matrix for speed
%         for bar=1:width:res %creating ternary white noise bars
%           if randseed(bar) <= p
%             barmat(:,bar:bar+width-1)=255;
%           elseif randseed(bar) >= 1-p
%             barmat(:,bar:bar+width-1)=0;
%           else
%             barmat(:,bar:bar+width-1)=gray;
%           end
%         end
%     end
% 
% seedmat.stims(i,:) = barmat(1,:);
% seedmat.stimdurs(i,1) = durseed;
% seedmat.p = p;
% seedmat.width = width;
% seedmat.tilt = tilt;
% 
% 
% outmat_temp=barmat; %only uses window of size of scaf
% 
% if i==1
%     outmat=outmat_temp;
% else
%     outmat=cat(3, outmat, outmat_temp);
% end
% 
% 
% end