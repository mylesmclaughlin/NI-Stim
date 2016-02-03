function M = NIStimMacro_tACS_SubthresAM


%% - Setup the marco struture -

% Setup the macro variables
M(1).macro.basefilename = 'RatTACS-'; % set base file name
M(1).macro.basefilename = [M(1).macro.basefilename  strrep(strrep(datestr(now),':','-'),' ','-')]; % add date to make filename unique
M(1).macro.nreps = 5; % - number of time to repeated each stimulus contained in the macro structure
M(1).macro.randomize = 1;  % randomize the order of stimului presentation
M(1).macro.record = 1; % record data (separate file for each rep of macro)

if M(1).macro.record == 0
    M(1).macro.allinonebinfile = 1; % specify if the recorded data will be stored in one continuous .bin file ...
elseif M(1).macro.record == 1
    M(1).macro.allinonebinfile = 0; %... or if each element in the array will have its own .bin file
end

% Set fixed variables for easy modification
numberreps = 1 ;
ampRange = [1:0.2:3];
freq = 300;
burstdur = 80;
burstrepperiod = 2000;
waveformindex = 1;
phase1pulsewidth = 300;
phase2pulsewidth = 300;

subthresAmp = 0.5;
subthresFreq = 2;
subthresCar = freq;
subthresAmpModDepth = 100;

%% --Setup stimulation variables in a struct array --

% set supra threshold same for all conditions
for n = 1:3
    M(n).stim.stim = 0;
    M(n).stim.amplitude = ampRange;
    M(n).stim.frequency = freq;
    M(n).stim.phase = 0;
    M(n).stim.continuous = 0;
    M(n).stim.ramp = 0;
    M(n).stim.burstdur = burstdur;
    M(n).stim.burstrepperiod = burstrepperiod;
    M(n).stim.numberreps = numberreps;
    M(n).stim.repsplayed = 0;
    M(n).stim.ampmoddepth = 0;
    M(n).stim.ampmodfreq = 2;
    M(n).stim.ampmodphase = 0;
    M(n).stim.waveformindex = waveformindex;
    M(n).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
    M(n).stim.phase1pulsewidth = phase1pulsewidth;
    M(n).stim.phase2pulsewidth = phase2pulsewidth;
    M(n).stim.phase1amp = 100;
    M(n).stim.phase2amp = -100;
    M(n).stim.phasegap = 0;
    M(n).stim.sameonallchannels = 1;
    M(n).stim.randomizesequence = 0;
end


% set continuous sub threshold different for each condition
n=0;
% Baseline stimulus
n=n+1;
M(n).stim.seqname = 'baseline'; %

% delay normal stimulus
M(n).stim.burstdelay = 0;
% settings for base stimulus
M(n).basestim.stim = 0;
M(n).basestim.amplitude = subthresAmp;
M(n).basestim.frequency = 2;
M(n).basestim.phase = 0;
M(n).basestim.ampmoddepth = 0;
M(n).basestim.ampmodfreq = 5;
M(n).basestim.ampmodphase = 0.5;
M(n).basestim.waveformindex = 1;
M(n).basestim.phase1pulsewidth = 50;
M(n).basestim.phase2pulsewidth = 50;
M(n).basestim.phase1amp = 100;
M(n).basestim.phase2amp = -100;
M(n).basestim.phasegap = 0;
M(n).basestim.waveformindex = 1;

% Single Frequency subthreshold stimulus
n=n+1;
M(n).stim.seqname = 'singlefreq'; %

% delay normal stimulus
M(n).stim.burstdelay = 1000*(1/subthresFreq)/4;
% settings for base stimulus
M(n).basestim.stim = 1;
M(n).basestim.amplitude = subthresAmp;
M(n).basestim.frequency = subthresFreq;
M(n).basestim.phase = 0;
M(n).basestim.ampmoddepth = 0;
M(n).basestim.ampmodfreq = 5;
M(n).basestim.ampmodphase = 0.5;
M(n).basestim.waveformindex = 1;
M(n).basestim.phase1pulsewidth = 50;
M(n).basestim.phase2pulsewidth = 50;
M(n).basestim.phase1amp = 100;
M(n).basestim.phase2amp = -100;
M(n).basestim.phasegap = 0;
M(n).basestim.waveformindex = 1;


% Amp Mod subthreshold stimulus
n=n+1;
M(n).stim.seqname = 'ampmod'; %

% delay normal stimulus
M(n).stim.burstdelay = 1000*(1/subthresFreq)/4 - 1000*(1/subthresCar)/4;% delay normal stimulus;
% settings for base stimulus
M(n).basestim.stim = 1;
M(n).basestim.amplitude = subthresAmp;
M(n).basestim.frequency = subthresCar;
M(n).basestim.phase = 0;
M(n).basestim.ampmoddepth = subthresAmpModDepth;
M(n).basestim.ampmodfreq = subthresFreq*2;
M(n).basestim.ampmodphase = 0.5;
M(n).basestim.waveformindex = 1;
M(n).basestim.phase1pulsewidth = 50;
M(n).basestim.phase2pulsewidth = 50;
M(n).basestim.phase1amp = 100;
M(n).basestim.phase2amp = -100;
M(n).basestim.phasegap = 0;
M(n).basestim.waveformindex = 1;



%% - Run the macro -
M = NIStimMacro(M);