function [ Hs2d ] = mk_hartleyset4dualstim( L, freqvec )
%mk_hartleyset4dualstim generate set of hartley stimuli (not rotated)
%   Detailed explanation goes here

% rotvec=0:30:330;
% freqvec=pi/6:pi/6:2*pi;

xs=(1:L)-L/2;
M=round(L/3);

for ff=1:length(freqvec)
   xs2=xs*pi*(freqvec(ff))/M;
   Hs(ff,:)=sin(xs2)+cos(xs2);
   Hs2d(ff,:,:)=repmat(Hs(ff,:),L,1);
end

end

