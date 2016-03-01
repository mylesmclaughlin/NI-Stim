function M = NIStimMacro_tDCS


%% - Setup the marco struture -

% Setup the macro variables
M(1).macro.basefilename = 'RatTDCSSubThres-'; % set base file name
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
ampRange = [1.8:.2:2.6];
freq = 300;
burstdur = 100;
burstrepperiod = 2000;
waveformindex = 2;
phase1pulsewidth = 200;
phase2pulsewidth = 200;

subthresAmp =[-2: .4 : 2];
subthresFreq = 32;
subthresCar = freq;
subthresAmpModDepth = 0;

%% --Setup stimulation variables in a struct array --

% set supra threshold same for all conditions
for n = 1:length(subthresAmp)
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

for n = 1:length(subthresAmp)
    M(n).stim.seqname = ['dc' num2str(subthresAmp(n))]; %
    
    % delay normal stimulus
    M(n).stim.burstdelay = 0; %1000*3*(1/subthresFreq)/4;
    % settings for base stimulus
    M(n).basestim.stim = 1;
    M(n).basestim.amplitude = subthresAmp(n);
    M(n).basestim.frequency = subthresFreq;
    M(n).basestim.phase = 0;
    M(n).basestim.dc = 1;
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
end



%% - Run the macro -
M = NIStimMacro(M);
