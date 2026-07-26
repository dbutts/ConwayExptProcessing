function [Kofiko_ET_TS, Kofiko_Xpix, Kofiko_Ypix, KofikoGains, KofikoOffsets, KofikoGains_Plexon, KofikoOffsets_Plexon] = loadKofikoEyeData(g_strcts,t_plexon, B)


%% Kofiko eye data
% Get the screen dimensions
ScreenSizeX_pix = g_strcts.g_strctStimulusServer.m_aiScreenSize(3);
ScreenSizeY_pix = g_strcts.g_strctStimulusServer.m_aiScreenSize(4);

% Convert Kofiko eye signal timestamps to Plexon time
Kofiko_ET_TS = g_strcts.g_strctEyeCalib.EyeRaw.TimeStamp';

% Get Kofiko raw eye signals
Kofiko_EyeRawX = g_strcts.g_strctEyeCalib.EyeRaw.Buffer(:,1);
Kofiko_EyeRawY = g_strcts.g_strctEyeCalib.EyeRaw.Buffer(:,2);
% Get Kofiko eye signal gains
Kofiko_GainX = g_strcts.g_strctEyeCalib.GainX.Buffer;
Kofiko_GainY = g_strcts.g_strctEyeCalib.GainY.Buffer;
% Get Kofiko eye signal gain timestamps & convert to plexon time
Kofiko_GainX_TS = g_strcts.g_strctEyeCalib.GainX.TimeStamp;

Kofiko_GainY_TS = g_strcts.g_strctEyeCalib.GainY.TimeStamp;

% Get Kofiko eye signal offsets
Kofiko_CenterX = g_strcts.g_strctEyeCalib.CenterX.Buffer;
Kofiko_CenterY = g_strcts.g_strctEyeCalib.CenterY.Buffer;
% Get Kofiko eye signal offset timestamps
Kofiko_CenterX_TS = g_strcts.g_strctEyeCalib.CenterX.TimeStamp;

Kofiko_CenterY_TS = g_strcts.g_strctEyeCalib.CenterY.TimeStamp;


% Determine the Kofiko gains and offsets for each eye signal time stamp

Kofiko_GainX_forEachTimeStamp = loadKofikoEyeData_helper(Kofiko_ET_TS, Kofiko_GainX_TS, Kofiko_GainX);
Kofiko_GainY_forEachTimeStamp = loadKofikoEyeData_helper(Kofiko_ET_TS, Kofiko_GainY_TS, Kofiko_GainY);
Kofiko_CenterX_forEachTimeStamp = loadKofikoEyeData_helper(Kofiko_ET_TS, Kofiko_CenterX_TS, Kofiko_CenterX);
Kofiko_CenterY_forEachTimeStamp = loadKofikoEyeData_helper(Kofiko_ET_TS, Kofiko_CenterY_TS, Kofiko_CenterY);


Kofiko_Xpix = Kofiko_GainX_forEachTimeStamp.*(Kofiko_EyeRawX - Kofiko_CenterX_forEachTimeStamp) + ScreenSizeX_pix/2;
Kofiko_Ypix = Kofiko_GainY_forEachTimeStamp.*(Kofiko_EyeRawY - Kofiko_CenterY_forEachTimeStamp) + ScreenSizeY_pix/2;

KofikoGains = [Kofiko_GainX_forEachTimeStamp, Kofiko_GainY_forEachTimeStamp];
KofikoOffsets = [Kofiko_CenterX_forEachTimeStamp, Kofiko_CenterY_forEachTimeStamp];

if nargin == 3

    Kofiko_GainX_TS_PlexonTime = [ones(size(Kofiko_GainX_TS')) Kofiko_GainX_TS']*B;
    Kofiko_GainY_TS_PlexonTime = [ones(size(Kofiko_GainY_TS')) Kofiko_GainY_TS']*B;
    Kofiko_CenterX_TS_PlexonTime = [ones(size(Kofiko_CenterX_TS')) Kofiko_CenterX_TS']*B;
    Kofiko_CenterY_TS_PlexonTime = [ones(size(Kofiko_CenterY_TS')) Kofiko_CenterY_TS']*B;

    Kofiko_GainX_forEachPlexonSample =loadKofikoEyeData_helper(t_plexon, Kofiko_GainX_TS_PlexonTime, Kofiko_GainX);
    Kofiko_GainY_forEachPlexonSample =loadKofikoEyeData_helper(t_plexon, Kofiko_GainY_TS_PlexonTime, Kofiko_GainY);
    Kofiko_CenterX_forEachPlexonSample =loadKofikoEyeData_helper(t_plexon, Kofiko_CenterX_TS_PlexonTime, Kofiko_CenterX);
    Kofiko_CenterY_forEachPlexonSample =loadKofikoEyeData_helper(t_plexon, Kofiko_CenterY_TS_PlexonTime, Kofiko_CenterY);

    KofikoGains_Plexon = [Kofiko_GainX_forEachPlexonSample',Kofiko_GainY_forEachPlexonSample'];
    KofikoOffsets_Plexon = [Kofiko_CenterX_forEachPlexonSample',Kofiko_CenterY_forEachPlexonSample'];

else
    KofikoGains_Plexon = [];
    KofikoOffsets_Plexon = [];
end





end


