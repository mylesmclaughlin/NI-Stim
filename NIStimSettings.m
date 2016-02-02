function S = NIStimSettings;

%% ----- NI settings ------
device = daq.getDevices;
if isempty(device)
    disp('No NI device found. Starting in playback mode');
    S.ni.connected = 0;
    S.ni.devname = '';
    S.ni.description = 'no device';
else
    if size(device,2) == 1
        devInd = 1;
    else
        devInd = [];
        for n = 1:length(device)
            if strcmp(device(n).Description,'National Instruments USB-6343')
                devInd = n;
                break
            end
        end
        if isempty(devInd)
            for n = 1:length(device)
                if strcmp(device(n).Description,'National Instruments USB-6216 (BNC)')
                    devInd = n;
                    break
                end
            end
        end
    end
    if isempty(devInd)
        disp('No NI device found. Starting in playback mode');
        S.ni.connected = 0;
        S.ni.devname = '';
        S.ni.description = 'no device';
    else
        disp(['Connected to ' device(devInd).Description])
        S.ni.connected = 1;
        S.ni.devname = device(devInd).ID;
        S.ni.description = device(devInd).Description;
    end
end

if strcmp(S.ni.description,'National Instruments USB-6216 (BNC)')
    disp(['Applying settings for ' S.ni.description] )
    S.ni.chin = [0 1 2 3 4 5];
    S.ni.chilabel = {'Current','Voltage','Trigger','X','Y','Z'};
    S.ni.chout = [0 1];
    S.ni.rate = 5e3;
    S.ni.voltrange =  [-1 1; -10 10; -10 10; -10 10; -10 10; -10 10;];
    S.ni.inputtype = {'SingleEnded','SingleEnded','SingleEnded','SingleEnded','SingleEnded','SingleEnded'};
elseif strcmp(S.ni.description,'National Instruments USB-6343')
    disp(['Applying settings for ' S.ni.description])
    S.ni.chin = [1 2]; %  [1 2 3 18 19 22];   ; %      %[1 2 3]; %  %[18 19 22 23]; % %  % % %
    S.ni.chilabel = {'Current','Voltage'}; %{'Current','Voltage','Trigger','X','Y','Z'}; %
    S.ni.chout = [0 1 2]; % 0 out is always trigger
    S.ni.rate = 200e3;
    S.ni.voltrange = [-1 1; -10 10];% [-1 1; -10 10; -10 10; -10 10; -10 10; -10 10;]; %[-0.5 0.5; -0.5 0.5; -0.5 0.5;]; %
    S.ni.inputtype = {'SingleEnded','SingleEnded'}; %{'SingleEnded','SingleEnded','SingleEnded','SingleEnded','SingleEnded','SingleEnded'}; %{'Differential','Differential','Differential'}; % %{'SingleEnded','Differential','SingleEnded','SingleEnded','SingleEnded','SingleEnded'}; %{'SingleEnded','Differential','SingleEnded'};
elseif strcmp(S.ni.description,'no device')
    S.ni.chin = [1 2 3 18 19 22];
    S.ni.chout = [0 1 2];
    S.ni.rate = 200e3;
    S.ni.voltrange = [-10 10];
    S.ni.chilabel = {'Data','Current','Voltage','Trigger'};
else % unknown device - try these settings
    disp(['No settings saved for ' S.ni.description '. Using generic NI settings'])
    S.ni.chin = [0 1];
    S.ni.chout = [0 1];
    S.ni.rate = 20e3;
    S.ni.voltrange = [-10 10];
end
S.ni.nchin = length(S.ni.chin);
S.ni.nchout = length(S.ni.chout);
S.ni.nupdatespersec =  2; % get data and update maximum of 5 times per second
S.ni.buffersize = S.ni.rate/S.ni.nupdatespersec;
S.ni.updateperiod = S.ni.buffersize/S.ni.rate;
if S.ni.updateperiod<0.1
    S = 'Update period is too fast';
    disp(S)
    return
end

%% ------ Add Paths ----
S.path = [pwd '\'];
addpath([pwd '\AccelCalib\'])
addpath([pwd '\BinaryFiles\'])
addpath([pwd '\General\'])
addpath([pwd '\Process\'])
addpath([pwd '\StimulusGeneration\'])
addpath([pwd '\Analysis\'])
addpath([pwd '\Macro\'])
addpath([pwd '\CurrentSource\'])
addpath([pwd '\CurrentSource\DS5\'])
addpath([pwd '\CurrentSource\AM2200\'])
% addpath(pwd)

%% ------ Amplifier Settings ------
S.amp.name = '';
S.amp.gain = 1;
S.amp.filter = [NaN NaN];

%% ------ Current Source Settings ------
S.current.present = 1;
S.current.namelist = {'DS5 10mA','DS5 2mA (Patients use 400kOhm cable)','AM 2200','None'};
S.current.namelistunits = {'mA','mA','mA','V'};
S.current.namevalue = 1;
S.current.onevoltequalsXmilliampslist = [1 0.2 0.1 1]; %0.1;
S.current.name = S.current.namelist(S.current.namevalue);
S.current.onevoltequalsXmilliamps = S.current.onevoltequalsXmilliampslist(S.current.namevalue);

%% ------ Data ------
S.data.dir = 'C:\Users\Public\Data\NIData\';
S.data.datatype = 'double'; %'single'
S.data.header = ['Fs=' num2str(S.ni.rate) ', Filter=[' num2str(S.amp.filter) '], Gain=' num2str(S.amp.gain)  ', Amp=' S.amp.name];

%% ------ Trigger -------
S.trigger.on = 1;
S.trigger.dur = 10; %ms
S.trigger.amp = 1;
S.trigger.chind = 3;

%% ------ Stimulate ------
S.stim.stim = 0;
S.stim.amplitude = 1;
S.stim.frequency = 100;
S.stim.phase = 0;
S.stim.continuous = 1;
S.stim.ramp = 1;
S.stim.burstdur = 100;
S.stim.burstrepperiod = 500;
S.stim.numberreps = Inf;
S.stim.repsplayed = 0;
S.stim.ampmoddepth = 0;
S.stim.ampmodfreq = 5;
S.stim.ampmodphase = 0.5;
S.stim.waveformindex = 1;
S.stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
S.stim.phase1pulsewidth = 50;
S.stim.phase2pulsewidth = 50;
S.stim.phase1amp = 100;
S.stim.phase2amp = -100;
S.stim.phasegap = 0;
S.stim.stimdir = 'C:\Users\u0043883\Google Drive\Work\MATLAB\NI-Stim\Stimuli\';
S.stim.customfilename = '';
S.stim.customdata = [];
S.stim.sameonallchannels = 1;
S.stim.stoppingmode = 0;
S.stim.randomizesequence = 0;


%% ----- Sequence Stimulus ----
S.sequence.on = 0;
S.sequence.thisseq = 0;
S.sequence.loopthisseq = S.sequence.thisseq;
S.sequence.nseq = 0;
S.sequence.parametername = '';
S.sequence.parametervalues = [];
S.sequence.data = [];
S.sequence.seq = [1:S.sequence.nseq];
S.sequence.seqIndex =  S.sequence.seq;

%% ----- Impedance Monitoring -----
S.impedance.calc = 0;
S.impedance.currentchind = 2;
S.impedance.voltagechind = 3;
S.impedance.currentconversion = 10e-3;
S.impedance.voltageconversion = 20; % DS5 use 20; SRS 560 Amp use Gain setting
S.impedance.val = NaN;

%% ------ Record ------
if S.ni.connected == 1
    S.rec.playbackmode = 0;
else
    S.rec.playbackmode = 1;
end
S.rec.rec = 0;
S.rec.plotdur = 5;
S.rec.maxplotfs = 10e3; % limit the amount of samples plotted to speed up system
if S.ni.rate <= S.rec.maxplotfs;
    S.rec.plotdownsample = 1;
    S.rec.plotfs = S.ni.rate;
else
    S.rec.plotdownsample = ceil(S.ni.rate/S.rec.maxplotfs);
    S.rec.plotfs = S.ni.rate/S.rec.plotdownsample;
end
S.rec.plotbuffersize = S.ni.buffersize/S.rec.plotdownsample;
S.rec.rawplotbuffersize = S.ni.buffersize;
S.rec.timevec = [-S.rec.plotdur:1/S.rec.plotfs:-1/S.rec.plotfs];
S.rec.rawtimevec = [-S.rec.plotdur:1/S.ni.rate:-1/S.ni.rate];
S.rec.plotdata = zeros(length(S.rec.timevec),S.ni.nchin);
S.rec.rawplotdata = zeros(length(S.rec.rawtimevec),S.ni.nchin);
S.rec.procplotdata = zeros(length(S.rec.rawtimevec),S.ni.nchin);
S.rec.ylims = [-10 10];
S.rec.xlims = [-S.rec.plotdur 0];
S.rec.filename = '';
S.rec.quickrec = 0;

S.rec.showplot = 1; % swithc off plotting to speed up system
S.rec.listen = 0;
S.rec.listenchind = 2;
S.rec.listendata = [];

%% ------ Accelerometer ------
S.accel.accel = 1;
if exist([S.data.dir 'AccelCalib\CalibrationData_' date '.mat'],'file')~=0
    load([S.data.dir 'AccelCalib\CalibrationData_' date '.mat']);
    S.accel.donecalib = 1;
    S.accel.callibration.g0 = g0;
    S.accel.callibration.g1 = g1;
else
    S.accel.donecalib = 0;
    S.accel.callibration.g0 = [0 0 0];
    S.accel.callibration.g1 = [0 0 0];
end
S.accel.bandpass = [0.25 50];
S.accel.filterorder = 2;
[S.accel.filterb S.accel.filtera] = butter(S.accel.filterorder,S.accel.bandpass/(S.ni.rate/2),'bandpass');
S.accel.chinds = [4 5 6];

%% ------ Process ------
S.proc.proc = 0;
S.proc.type = 'accel'; %'trigAvg';  %'trigSpike'; %   %
eval(['S = NIprocess_' S.proc.type '(''settings'',S);']);

%% ------ Macro ------
S.macro.active = 0;