function M = NIStimMacro_tACS_AMDepth


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
numberreps = 10;
ampRange = 1.5;
freqvec = 400;
burstdur = 2000;
burstrepperiod = 2000;
ampmoddepth = [0:20:100];
ampmodfreq = 4;
waveformindex = 1; % 1 for sin 2 for pulse
phase1pulsewidth = 200;
phase2pulsewidth = 200;

% --Setup stimulation variables in a struct array --
% Stimulus 1 - Sin 50 Hz
for n = 1:length(ampmoddepth) %
    M(n).stim.seqname = [num2str(ampmoddepth(n))];
    M(n).stim.stim = 0;
    M(n).stim.amplitude = ampRange;
    M(n).stim.frequency = freqvec;
    M(n).stim.phase = 0;
    M(n).stim.continuous = 0;
    M(n).stim.ramp = 0;
    M(n).stim.burstdur = burstdur;
    M(n).stim.burstrepperiod = burstrepperiod;
    M(n).stim.numberreps = numberreps;
    M(n).stim.repsplayed = 0;
    M(n).stim.ampmoddepth = ampmoddepth(n);
    M(n).stim.ampmodfreq = ampmodfreq;
    M(n).stim.ampmodphase = 0;
    M(n).stim.waveformindex = waveformindex;
    M(n).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
    M(n).stim.phase1pulsewidth = phase1pulsewidth;
    M(n).stim.phase2pulsewidth = phase2pulsewidth;
    M(n).stim.phase1amp = 100;
    M(n).stim.phase2amp = -100;
    M(n).stim.phasegap = 0;
    M(n).stim.sameonallchannels = 1;
end

%% - Run the macro -
M = NIStimMacro(M);
