Robs = data(1).Robs;
stim = data(1).stim;

nSpks = sum(Robs,2);
nFrames = size(Robs,2);

nLags = 10;

%% by time
propSpikes = [1, 99/100, 19/20, 9/10, 3/4, 2/3, 1/2, 1/3, 1/4, 1/10, 1/20, 1/100];
for i = 1:numel(propSpikes)
    nSampleFrames = round(propSpikes(i) .* nFrames); %round(interp1(x_sub, y_sub, propSpikes(i), 'spline'));
    % idx = sort(randperm(nFrames, nSampleFrames));
    idx = 1:nSampleFrames;
    STA_sub{i} = generate_stas(Robs(:,idx),stim(:,:,:,idx), 10);
    truePropSpikes(:,i) = sum(Robs(:,idx),2)./nSpks;
    corrMat(:,:,i) =  corr(transpose(reshape(STA_sub{1}.DKL, 43, [])),  transpose(reshape(STA_sub{i}.DKL, 43, [])));
    c(:,i) = diag(corrMat(:,:,i));
    sprintf('i/ %i\n', numel(propSpikes));

end

figure, hold on
subplot(1,3,1)
plot(truePropSpikes', c', 'ko-');
ylabel('Correlation')
xlabel('Proportion of spikes')

subplot(1,3,2)
plot(truePropSpikes' .* nSpks', c', 'ko-')
xlabel('Number of spikes')

subplot(1,3,3)
fps = 60;
ts = round(propSpikes .* nFrames) / fps;
plot(ts', c', 'ko-')
xlabel('Time (s)')




%% discontinous/random sampling

propSpikes = [1, 99/100, 19/20, 9/10, 3/4, 2/3, 1/2, 1/3, 1/4, 1/10, 1/20, 1/100];
for i = 1:numel(propSpikes)
    nSampleFrames = round(propSpikes(i) .* nFrames); %round(interp1(x_sub, y_sub, propSpikes(i), 'spline'));
    idx = sort(randperm(nFrames, nSampleFrames));
    %idx = 1:nSampleFrames;
    STA_sub2{i} = generate_stas(Robs(:,idx),stim(:,:,:,idx), 10);
    truePropSpikes2(:,i) = sum(Robs(:,idx),2)./nSpks;
    corrMat2(:,:,i) =  corr(transpose(reshape(STA_sub2{1}.DKL, 43, [])),  transpose(reshape(STA_sub2{i}.DKL, 43, [])));
    c2(:,i) = diag(corrMat2(:,:,i));
    sprintf('i/ %i\n', numel(propSpikes));

end

figure, hold on
subplot(1,3,1)
plot(truePropSpikes2', c2', 'ko-');
ylabel('Correlation')
xlabel('Proportion of spikes')

subplot(1,3,2)
plot(truePropSpikes2' .* nSpks', c2', 'ko-')
xlabel('Number of spikes')

subplot(1,3,3)
fps = 60;
ts = round(propSpikes .* nFrames) / fps;
plot(ts', c2', 'ko-')
xlabel('Time (s)')

for i = 1:10
    nSampleFrames = round(0.8 .* nFrames);
    idx = sort(randperm(nFrames, nSampleFrames));
    STA80{i} = generate_stas(Robs(:,idx),stim(:,:,:,idx), 10);
end


%%

stimMat = permute(data(1).stim, [4 1 2 3]);

stimMat = reshape(stimMat, size(data(1).stim, 4), []);

[C, ia, ic] = unique(stimMat, 'rows');



for lag = 0:9
    for i = 1:size(C,1)
        idx = find(ic == i) + lag;
        idx(idx > size(Robs,2)) = [];
        numSpks(:,i, lag+1) = sum(Robs(:,idx), 2);
        meanSpks(:,i,lag+1) =  mean(Robs(:,idx), 2);
        stdSpks(:,i,lag+1) = std(Robs(:,idx),[], 2);
    end
end

lag = 4;
unit = 19;
n = size(C,1);
[B,I] = sort(numSpks(unit,:, lag+1), 'descend');
nIms = single(C(I(1:n),:));
W = B(1:n);
recon = reshape((mean(W'.*nIms,1)), 60,60,3);
recon = circshift(recon,[30,30]);
figure, imagesc(recon(:,:,1)); colormap gray

idx = (size(C,1)-n+1):size(C,1);
W_inv = 1./(eps + B(idx));
nIms_inv = single(C(I(idx),:));
recon_inv = reshape((mean(W_inv'.*nIms_inv,1)), 60,60,3);

recon_inv= circshift(recon_inv,[30,30]);
figure, imagesc(recon_inv(:,:,1)); colormap gray

% dumb eye position correction?

h = squeeze(STA.DKL(unit,:,:,1,lag+1));
for i = 1:size(C,1)
    im = reshape(C(I(i),1:3600), 60,60);
    fim = imfilter(im,h);
    [atemp, btemp] = find(fim == max(fim(:)));
    rr(i) = atemp(1); cc(i) = btemp(1);
end
for i = 1:size(C,1)
    temp = reshape(C(I(i),:), 60,60,3);
    temp = circshift(temp, [-rr(i), -cc(i)]);
    temp = reshape(temp, 10800, []);
    C_shifted(i,:) = temp;
end
nIms = single(C_shifted(I(1:n),:));
W = B(1:n);
recon = reshape((mean(W'.*nIms,1)), 60,60,3);
%recon = circshift(recon,[30,30]);
figure, imagesc(recon(:,:,1)); colormap gray