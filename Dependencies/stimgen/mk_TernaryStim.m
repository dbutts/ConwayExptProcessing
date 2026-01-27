function [outmat, seedmat]=mk_TernaryStim(p, width, tilt, stim_frames, res, seedmat_in)
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

if p>.5; error('p is too large'); end
if nargin < 3; tilt=0; end
if nargin < 4; stim_frames=1; end
if nargin < 5; res=80; end

if nargin == 6 % Loads in information from seed matrix
p = seedmat_in.p;
width = seedmat_in.width;
tilt = seedmat_in.tilt;
stim_size = seedmat_in.stim_size;
end

iterations=stim_frames; %ceil(run_duration*60/mean(stim_duration)); %turns run duration from seconds into frames

for i=1:iterations %looping through iterations; stimulus duration is added

    if nargin == 6
        durseed = seedmat_in.stimdurs(i,1);
        barmat = repmat(seedmat_in.stims(i,:), res, 1);
    else
        randseed=rand(1,res);

        if length(stim_duration)==1
            durseed=stim_duration;
        elseif length(stim_duration)==2
            mindur=stim_duration(1);
            maxdur=stim_duration(2);
            durseed=randi([mindur maxdur]); %generates random spread for durations
        else
            error('Stimulus Duration input must be 1 or 2 inputs')
        end

        barmat=zeros(res, res); %preallocate matrix for speed
        for bar=1:width:res %creating ternary white noise bars
          if randseed(bar) <= p
            barmat(:,bar:bar+width-1)=255;
          elseif randseed(bar) >= 1-p
            barmat(:,bar:bar+width-1)=0;
          else
            barmat(:,bar:bar+width-1)=gray;
          end
        end
    end

seedmat.stims(i,:) = barmat(1,:);
seedmat.stimdurs(i,1) = durseed;
seedmat.p = p;
seedmat.width = width;
seedmat.tilt = tilt;


outmat_temp=barmat; %only uses window of size of scaf

if i==1
    outmat=outmat_temp;
else
    outmat=cat(3, outmat, outmat_temp);
end


end