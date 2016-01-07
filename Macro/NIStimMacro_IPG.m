function M = NIStimMacro_IPG


%% - Setup the marco struture -

% Setup the macro variables
n=1;
M(n).macro.basefilename = 'Rat-PMIPG-'; % set base file name
M(n).macro.basefilename = [M(n).macro.basefilename  strrep(strrep(datestr(now),':','-'),' ','-')]; % add date to make filename unique
M(n).macro.nreps = 5; % - number of time to repeated each stimulus contained in the macro structure
M(n).macro.randomize = 1;  % randomize the order of stimului presentation
M(n).macro.record = 0; % record data (separate file for each rep of macro)

if M(n).macro.record == 0
    M(n).macro.allinonebinfile = 1; % specify if the recorded data will be stored in one continuous .bin file ...
elseif M(n).macro.record == 1
    M(n).macro.allinonebinfile = 0; %... or if each element in the array will have its own .bin file
end

% Set fixed variables for easy modification
ampRange = 0.8:0.1:1.8;
frequency = 300;
burstdur = 100;
burstrepperiod = 2000;
phase1pulsewidth = 100;
phase2pulsewidth = 1000;
phase1amp = 100;
phase2amp = -10;

phasegap = [0 50 100 200 400];

% --Setup stimulation variables in a struct array --
% Stimulus 1 - IPG
for n = 1:length(phasegap)
    M(n).stim.seqname = ['IPG' num2str(phasegap(n))];
    M(n).stim.stim = 0;
    M(n).stim.amplitude = ampRange;
    M(n).stim.frequency = frequency;
    M(n).stim.phase = 0;
    M(n).stim.continuous = 0;
    M(n).stim.ramp = 0;
    M(n).stim.burstdur = burstdur;
    M(n).stim.burstrepperiod = burstrepperiod;
    M(n).stim.numberreps = 1;
    M(n).stim.repsplayed = 0;
    M(n).stim.ampmoddepth = 0;
    M(n).stim.ampmodfreq = 5;
    M(n).stim.ampmodphase = 0;
    M(n).stim.waveformindex = 2;
    M(n).stim.waveformlist = {'sine','pulse','triangle','gaussian','custom'}; % pulse, custom
    M(n).stim.phase1pulsewidth = phase1pulsewidth;
    M(n).stim.phase2pulsewidth = phase2pulsewidth;
    M(n).stim.phase1amp = phase1amp;
    M(n).stim.phase2amp = phase2amp;
    M(n).stim.phasegap = phasegap(n);
    M(n).stim.sameonallchannels = 1;
    M(n).stim.randomizesequence = 0;
end

%% - Run the macro -
M = NIStimMacro(M);
