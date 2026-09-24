function PlexET_ad_calib = loadPlexonEyeData(plexon_fname, chanNums, rig, plexonAnalogScale, gains, ET_Eyelink)


if nargin < 6
    ET_Eyelink = 1;
    bino_ddpi = 1;
else
    bino_ddpi = 0;
end

if nargin < 5 % if gains not provided
    gains = [1 1];
end

if nargin < 4 % if plexon scale not provided
    plexonAnalogScale = 1;
end

pl2 = PL2ReadFileIndex(plexon_fname);
temp = vertcat(pl2.AnalogChannels{:});
analogChanNames = {temp.Name};
numDigitsInLastAIchan = ceil(log10(sum(contains(analogChanNames, 'AI'))));

for i = 1:numel(chanNums)
    chanName =  ['AI' num2str(chanNums(i), ['%0' num2str(numDigitsInLastAIchan) '.f'])];
    [~,~,~,~,AI{i}] = plx_ad_v(plexon_fname, chanName);
end

sync_ad = AI{1};
arc_ad = AI{2};


if strcmpi(rig, 'B')
    if bino_ddpi

        eyeX1_ad = AI{3}; % left
        eyeY1_ad = AI{4}; % left
        eyeX2_ad = AI{5}; % right
        eyeY2_ad = AI{6}; % right
        pupil1_ad = AI{7}; % left
        pupil2_ad = AI{8}; % right

        % calibrate

        eyeX1_ad_calib = plexonAnalogScale .* gains(:,1).* eyeX1_ad';
        eyeY1_ad_calib = plexonAnalogScale .* gains(:,2).* eyeY1_ad';

        eyeX2_ad_calib = plexonAnalogScale .* gains(:,1).* eyeX2_ad';
        eyeY2_ad_calib = plexonAnalogScale .* gains(:,2).* eyeY2_ad';

    else
    end
elseif strcmpi(rig, 'C')
    if bino_ddpi

        eyeX2_ad = AI{5}; % right
        eyeY2_ad = AI{6}; % right
        eyeX1_ad = AI{7}; % left
        eyeY1_ad = AI{8}; %left
        pupil1_ad = AI{3}; % left
        pupil2_ad = AI{4}; % right

        eyeX1_ad_calib = plexonAnalogScale .* gains(:,1) .*eyeX1_ad;
        eyeY1_ad_calib = plexonAnalogScale .* gains(:,2) .*eyeY1_ad;

        eyeX2_ad_calib = plexonAnalogScale .* gains(:,1) .*eyeX2_ad;
        eyeY2_ad_calib = plexonAnalogScale .* gains(:,2) .*eyeY2_ad;

    elseif ET_Eyelink == 1
        eyeX1_ad = AI{5};
        eyeY1_ad = AI{6};
        eyeX2_ad = AI{7};
        eyeY2_ad = AI{8};

        eyeX1_ad_calib = (eyeX1_ad' - plexonAnalogScale.*median(eyeX1_ad)) .* gains(:,1);
        eyeY1_ad_calib = (eyeY1_ad' - plexonAnalogScale.*median(eyeY1_ad)) .* gains(:,2);

        eyeX2_ad_calib = (eyeX2_ad' - plexonAnalogScale.*median(eyeX2_ad)) .* gains(:,1);
        eyeY2_ad_calib = (eyeY2_ad' - plexonAnalogScale.*median(eyeY2_ad)) .* gains(:,2);

        pupil1_ad = nan(size(eyeX1_ad_calib));
        pupil2_ad = nan(size(eyeX1_ad_calib));

    elseif ET_Eyelink == 2
        eyeX1_ad = AI{7};
        eyeY1_ad = AI{8};

        eyeX1_ad_calib = (eyeX1_ad' - plexonAnalogScale.*median(eyeX1_ad)) .* gains(:,1);
        eyeY1_ad_calib = (eyeY1_ad' - plexonAnalogScale.*median(eyeY1_ad)) .* gains(:,2);

        eyeX2_ad_calib = nan(size(eyeX1_ad_calib));
        eyeY2_ad_calib = nan(size(eyeX1_ad_calib));

        pupil1_ad = nan(size(eyeX1_ad_calib));
        pupil2_ad = nan(size(eyeX1_ad_calib));
    elseif ET_Eyelink == 3
        eyeX1_ad = AI{7};
        eyeY1_ad = AI{8};

        eyeX1_ad_calib = plexonAnalogScale .* gains(:,1) .*eyeX1_ad;
        eyeY1_ad_calib = plexonAnalogScale .* gains(:,2) .*eyeY1_ad;

        eyeX2_ad_calib = nan(size(eyeX1_ad_calib));
        eyeY2_ad_calib = nan(size(eyeX1_ad_calib));

        pupil1_ad = nan(size(eyeX1_ad_calib));
        pupil2_ad = nan(size(eyeX1_ad_calib));

    elseif ET_Eyelink == 0
        eyeX1_ad = AI{7};
        eyeY1_ad = AI{8};

        eyeX1_ad_calib = eyeX1_ad - median(eyeX1_ad);
        eyeY1_ad_calib = eyeY1_ad - median(eyeY1_ad);


        eyeX2_ad_calib = nan(size(eyeX1_ad_calib));
        eyeY2_ad_calib = nan(size(eyeX1_ad_calib));

        pupil1_ad = nan(size(eyeX1_ad_calib));
        pupil2_ad = nan(size(eyeX1_ad_calib));
    else
    end
end

PlexET_ad_calib(:,1) = sync_ad;
PlexET_ad_calib(:,2) = arc_ad;
PlexET_ad_calib(:,3) = pupil1_ad; % left
PlexET_ad_calib(:,4) = pupil2_ad; % right
PlexET_ad_calib(:,5) = eyeX2_ad_calib; % right
PlexET_ad_calib(:,6) = eyeY2_ad_calib; % right
PlexET_ad_calib(:,7) = eyeX1_ad_calib; % left
PlexET_ad_calib(:,8) = eyeY1_ad_calib; % left

end