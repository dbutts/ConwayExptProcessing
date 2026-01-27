function [RFarea,RFsize, fullRF] = fn_getRFarea(cur_filt_map,ratio)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    thresh=max(cur_filt_map(:))*ratio;
    B{1}=0;
    
    while length(find(cur_filt_map>thresh))>numel(cur_filt_map)*0.5
        ratio=ratio+0.1;
        thresh=max(cur_filt_map(:))*ratio;
    end

    cur_filt_bin = ones(size(cur_filt_map));
    cur_filt_bin(cur_filt_map<thresh)=0;

    B = bwboundaries(cur_filt_bin); rfbound=[];
    for iii=1:length(B); rfbound(iii)=length(B{iii});end
        [~,peakRF]=max(rfbound);
    
    RFarea = polyarea(B{peakRF}(:,1),B{peakRF}(:,2) );
    RFsize = sqrt(RFarea);
    fullRF = B;

end