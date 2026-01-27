%%
%samples = adReadWrapper(fullFileName, spkcOnChanNames{i+opts.ChnOffset})';

fidout = fopen(datPath, 'a'); % opening file for appending data
fwrite(fidout, samples, 'int16'); % append data
fclose(fidout); % close for memory

%%
datPath='/home/felix/Dropbox/Project_BevilColor/Simulations/KSdebug/appended1.dat'
samples=[1:20]
fidout = fopen(datPath, 'a'); % opening file for appending data
fwrite(fidout, samples, 'int16'); % append data
fclose(fidout); % close for memory

samples=[21:40]
fidout = fopen(datPath, 'a'); % opening file for appending data
fwrite(fidout, samples, 'int16'); % append data
fclose(fidout); % close for memory

fidout = fopen(datPath); % opening file for appending data
test1= fread(fidout)
%%
datPath='/home/felix/Dropbox/Project_BevilColor/Simulations/KSdebug/appended2.dat'
samples=[1:20]'
fidout = fopen(datPath, 'a'); % opening file for appending data
fwrite(fidout, samples, 'int16'); % append data
fclose(fidout); % close for memory

samples=[21:40]'
fidout = fopen(datPath, 'a'); % opening file for appending data
fwrite(fidout, samples, 'int16'); % append data
fclose(fidout); % close for memory

fidout = fopen(datPath); % opening file for appending data
test2= fread(fidout)

%%
datPath='/home/felix/Dropbox/Project_BevilColor/Simulations/KSdebug/appended0.dat'
samples=[1:20; 21:40; 41:60; 61:80]
fidout = fopen(datPath, 'w'); % opening file for appending data
fwrite(fidout, samples, 'int16'); % append data
fclose(fidout); % close for memory

fidout = fopen(datPath); % opening file for appending data
test3= fread(fidout)

%%
datPath='/home/felix/Dropbox/Project_BevilColor/Simulations/KSdebug/appendedtest1.dat'
samples=[1:20]';
fidout = fopen(datPath, 'w'); % opening file for appending data
fwrite(fidout, samples, 'int16', 6); % append data
fclose(fidout); % close for memory

fidout = fopen(datPath, 'r+'); % opening file for appending data
fseek(fidout, 2,"bof")
samples=[21:40]';
fwrite(fidout, samples, 'int16', 6); % append data
fclose(fidout); % close for memory

fidout = fopen(datPath, 'r+'); % opening file for appending data
fseek(fidout, 4,"bof")
samples=[41:60]';
fwrite(fidout, samples, 'int16', 6); % append data
fclose(fidout); % close for memory

fidout = fopen(datPath, 'r+'); % opening file for appending data
fseek(fidout, 6,"bof")
samples=[61:80]';
fwrite(fidout, samples, 'int16', 6); % append data
fclose(fidout); % close for memory

fidout = fopen(datPath); % opening file for appending data
test4= fread(fidout)

%%
datPath2='/home/felix/Dropbox/Project_BevilColor/Simulations/KSdebug/appendedtestnew.dat'
fidout2 = fopen(datPath2, 'w');
fidout = fopen(datPath, 'r+'); % opening file for appending data

fseek(fidout, 6, 'bof')

A = fread(fidout); %// Read the data
fwrite(fidout2, A); %// Write data to new file

%// Close the files
fclose(fidout);
fclose(fidout2);

fidout = fopen(datPath2); % opening file for appending data
test5= fread(fidout)
