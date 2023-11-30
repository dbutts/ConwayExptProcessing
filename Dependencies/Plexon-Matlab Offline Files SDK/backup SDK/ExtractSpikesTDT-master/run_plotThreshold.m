figure;
hs(1) = subplot(2,1,1);
plot(fdata(1,:));
hold on;
plot(ts,fdata(1,ts),'o');
hs(2) = subplot(2,1,2);
plot(SNLEdata(1,:));
hold on;
plot([0 length(SNLEdata(1,:))],[thresholds(1) thresholds(1)],'k-');

linkaxes(hs,'x');