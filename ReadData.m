% This function reads in the data files and extracts the raw acceleration (in
% Volts), sampling frequency and time. 

% This function needs the full path of the datafile as input and gives the
% structure S as an output.S contains the raw acceleration
% (S.dataAccelRaw), sampling frequency (S.fs) and time vector (S.time)

function  S = ReadData(fileName)

% settings
accelChInd = [4 5 6];

% load data
D = binread(fileName);
commaInd = strfind(D.header,',');
equalInd = strfind(D.header,'=');
fs = str2num(D.header(equalInd(1)+1:commaInd(1)-1));
dataAccel = D.data;
dataAccelRaw = dataAccel(:,accelChInd);
current = dataAccel(:,1);

S.dataAccelRaw = dataAccelRaw;
S.fs = fs;
S.time  = [1:length(S.dataAccelRaw)]/S.fs;
S.current = current;