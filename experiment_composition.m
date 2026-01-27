function experiment_composition( exptdata )
% 
% Usage: experiment_composition( exptdata )
%
% Displays different types of stimuli used across all trials of given data

fprintf( '%d trials total\n', length(exptdata))

Ntr = length(exptdata);
stims = zeros(Ntr,2);

for tt=1:length(exptdata)
	stims(tt,1) = exptdata{tt, 1}.DualstimPrimaryuseRGBCloud; 
	stims(tt,2) = exptdata{tt, 1}.DualstimSecondaryUseCloud; 
	%test(tt)=exptdata{tt, 1}.m_bMonkeyFixated;
	%test(tt)=exptdata{tt, 1}.m_aiStimulusArea; 
	%test(tt)=exptdata{tt, 1}.m_aiStimulusRect(1);
	%test(tt)=exptdata{tt, 1}.m_aiStimulusRect(2);
end

types = unique(stims(:,1));
fprintf('Laminar probe stimuli types: ')
for ii=1:length(types)
	fprintf('%d ', types(ii))
end
fprintf('\n')

for ii=1:length(types)
	fprintf( '  %d: %d', types(ii), sum(stims(:,1) == types(ii)))
	switch(types(ii))
		case 8, disp(' (color cloud)');
		case 7, disp(' (luminance cloud)');
		case 6, disp(' (color Hartleys)')
		case 0, disp(' (ground truth)')
		otherwise 
			disp(' (Not sure)')
	end
end

types = unique(stims(:,2));
fprintf('\n\nET stimuli types: ')
for ii=1:length(types)
	fprintf('%d ', types(ii))
end
fprintf('\n')
for ii=1:length(types)
	fprintf( '  %d: %d', types(ii), sum(stims(:,2) == types(ii)))
	switch(types(ii))
		case 7, disp(' (color cloud)');
		case 6, disp(' (luminance cloud)');
		case 1, disp(' (1-D H-V alternating)')
		otherwise 
			disp(' (Not sure)')
	end
end
fprintf('\n')
