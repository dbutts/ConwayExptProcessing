function bufferVal_forEachTimePoint= loadKofikoEyeData_helper(all_ts, event_ts, bufferVals)

bufferVal_forEachTimePoint = nan(size(all_ts));

for i = 1:numel(bufferVals)
    idx = all_ts >= event_ts(i);
    bufferVal_forEachTimePoint(idx) = bufferVals(i);

end

end