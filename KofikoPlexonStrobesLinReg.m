function B = KofikoPlexonStrobesLinReg(g_strcts, plexon_fname)

% Read kofiko strobes from plexon and kofiko files
[events.count, events.timeStamps, events.strobeNumber] = plx_event_ts(plexon_fname, 257);

kofikoSyncStrobesTS = transpose(g_strcts.g_strctDAQParams.LastStrobe.TimeStamp(g_strcts.g_strctDAQParams.LastStrobe.Buffer == g_strcts.g_strctSystemCodes.m_iSync));
plexonSyncStrobesTS = events.timeStamps(events.strobeNumber == g_strcts.g_strctSystemCodes.m_iSync);

% Linear regression to get kofiko time stamps into plexon time
B = [ones(size(kofikoSyncStrobesTS)) kofikoSyncStrobesTS]\plexonSyncStrobesTS;

end