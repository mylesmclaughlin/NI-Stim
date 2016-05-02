function NIStimPreset(command)

if isstruct(command)
    P = command;
    command = 'init';
end

switch command
    case 'init'
        NISPinit(P);
    case 'stop'
        NISPstop
end

%--------------------------------------------------------------------------
function NISPinit(P);
global S

NISPupdateGUI(P)
if P.preset.record == 1
    NIStim('stimRec'); % record and stimulate
elseif P.preset.record == 0
    NIStim('startStim'); % stimulate only
end
disp(['Playing stimulus ' P.preset.basefilename])

% save preset datafile
save([S.data.dir P.preset.basefilename '.NISmacro'],'P','-mat');

while S.sequence.thisseq < length(S.sequence.seqIndex)
    pause(1)
end
pause(2*P.stim.burstrepperiod/1e3 + 3*S.ni.updateperiod) % wait for rep to end and update to happen before stopping
NIStim('stopStim')
NIStim('stopRec')

disp(['Finished stimulus ' P.preset.basefilename]) 

% %--------------------------------------------------------------------------
% function NISPstop
% global S
% S.macro.active = 0;

%--------------------------------------------------------------------------
function NISPupdateGUI(P)
% update Stim GUI
global S
set(S.stim.ampbut,'string',num2str(P.stim.amplitude));
set(S.stim.freqbut,'string',num2str(P.stim.frequency));
set(S.stim.phasebut,'string',num2str(P.stim.phase'));
set(S.stim.contbut,'value',P.stim.continuous);
set(S.stim.rampbut,'value',P.stim.ramp);
set(S.stim.burstdurbut,'string',P.stim.burstdur);
set(S.stim.burstrepperiodbut,'string',num2str(P.stim.burstrepperiod));
set(S.stim.numberrepsbut,'string',num2str(P.stim.numberreps));
set(S.stim.ampmoddepthbut,'string',num2str(P.stim.ampmoddepth));
set(S.stim.ampmodfreqbut,'string',num2str(P.stim.ampmodfreq));
set(S.stim.samechanbut,'value',P.stim.sameonallchannels);
set(S.stim.waveformbut,'value',P.stim.waveformindex);
set(S.stim.phase1pwbut,'string',num2str(P.stim.phase1pulsewidth));
set(S.stim.phase2pwbut,'string',num2str(P.stim.phase2pulsewidth));
set(S.stim.phase1ampbut,'string',num2str(P.stim.phase1amp));
set(S.stim.phase2ampbut,'string',num2str(P.stim.phase2amp));
set(S.stim.phasegapbut,'string',num2str(P.stim.phasegap));
S.stim.randomizesequence = P.stim.randomizesequence;

if isfield(P.stim,'series')
    S.stim.series.amplitude = P.stim.series.amplitude;
    S.stim.series.burstrepperiod = P.stim.series.burstrepperiod;
end

if isfield(P,'basestim')
      % delay normal stimulus
        S.stim.burstdelay = P.stim.burstdelay;
        % settings for base stimulus
        S.basestim.stim = P.basestim.stim;
        S.basestim.amplitude = P.basestim.amplitude;
        S.basestim.frequency = P.basestim.frequency;
        S.basestim.phase = P.basestim.phase;
        S.basestim.burstdur = P.basestim.burstdur;
        S.basestim.dc = P.basestim.dc;
        S.basestim.ampmoddepth = P.basestim.ampmoddepth;
        S.basestim.ampmodfreq = P.basestim.ampmodfreq;
        S.basestim.ampmodphase = P.basestim.ampmodphase;
        S.basestim.waveformindex = P.basestim.waveformindex;
        S.basestim.phase1pulsewidth = P.basestim.phase1pulsewidth;
        S.basestim.phase2pulsewidth = P.basestim.phase2pulsewidth;
        S.basestim.phase1amp = P.basestim.phase1amp;
        S.basestim.phase2amp = P.basestim.phase2amp;
        S.basestim.phasegap = P.basestim.phasegap;
        S.basestim.waveformindex = P.basestim.waveformindex;
end

set(S.rec.filenamebut,'string',P.preset.basefilename);
set(S.stim.filenamebut,'string',P.preset.basefilename);

% NIStim('parseGUI')
%NIStim('continuous')
%NIStim('waveform')
%NIStim('pulsewidth')
%NIStim('frequency')
