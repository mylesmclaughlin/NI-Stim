function M = NIStimMacro(command)

if isstruct(command)
    M = command;
    command = 'init';
end

switch command
    case 'init'
        M = NISMinit(M);
    case 'stop'
        NISMstop
    case 'pause'
        NISMpause
    case 'resume'
        NISMresume
    case 'playnext'
        NISMplaynext
end

%--------------------------------------------------------------------------
function M = NISMinit(M)
global S 

% Make randomzed macro struct
o=0;
if M(1).macro.randomize == 1
    for n = 1:M(1).macro.nreps
        rOrder = randperm(length(M));
        for m = 1:length(rOrder)
            o = o+1;
            N(o) = M(rOrder(m));
            N(o).macro = M(1).macro;
            N(o).macro.sequence = rOrder(m);
            N(o).macro.rep = n;
            N(o).macro.seqname = N(o).stim.seqname;
            N(o).macro.localfilename = [M(1).macro.basefilename '-' N(o).macro.seqname '-Seq' num2str(N(o).macro.sequence) '-Rep' num2str(n)];
            N(o).macro.done = 0;
        end
    end
end
M = N;
S.macro.basefilename = M(1).macro.basefilename;

% make macro fig
NISMmakeInputBox;

% save macro file
save([S.data.dir S.macro.basefilename '.NISmacro'],'M','-mat')

% run the macro
S.macro.run = 1; 
S.macro.active = 1; 
disp('Running macro...')

for n = 1:length(M)
    NISMplaynext;
    while S.sequence.thisseq < length(S.sequence.seqIndex)
        pause(1)
    end
    pause(2*M(n).stim.burstrepperiod/1e3 + 3*S.ni.updateperiod) % wait for rep to end and update to happen before stopping
    NIStim('stopStim')
    NIStim('stopRec')

    
    disp(['Recorded data to ' M(n).macro.localfilename])
    disp('--------------------------------------------')

    % check for pause
    while S.macro.run == 0
        pause(0.5)
        if S.macro.active == 0
            break
        end
    end
    
    % check for stop
    if S.macro.active == 0
        break
    end
end
disp('Marco has finished')

%--------------------------------------------------------------------------
function NISMplaynext
global S

% load macro data file
load([S.data.dir S.macro.basefilename '.NISmacro'],'M','-mat');
for i = 1:length(M)
    done(i) = M(i).macro.done;
end
zind = find(done==0);
n = zind(1);

% update NIStim gui
NISMupdateGUI(M(n));

% give system time to recover
pause(1)

% Stim & rec
if M(n).macro.record == 1
    NIStim('stimRec'); % record and stimulate
elseif M(n).macro.record == 0
    NIStim('startStim'); % stimulate only
end
disp(['Playing stimulus ' M(n).macro.localfilename])

% marks as done
M(n).macro.done = 1;

% update macro data file
save([S.data.dir M(1).macro.basefilename '.NISmacro'],'M','-mat');

if n == length(M)
    NISMstop
end

%--------------------------------------------------------------------------
function NISMupdateGUI(M)
% update Stim GUI
global S
set(S.stim.ampbut,'string',num2str(M.stim.amplitude));
set(S.stim.freqbut,'string',num2str(M.stim.frequency));
set(S.stim.phasebut,'string',num2str(M.stim.phase'));
set(S.stim.contbut,'value',M.stim.continuous);
set(S.stim.rampbut,'value',M.stim.ramp);
set(S.stim.burstdurbut,'string',M.stim.burstdur);
set(S.stim.burstrepperiodbut,'string',num2str(M.stim.burstrepperiod));
set(S.stim.numberrepsbut,'string',num2str(M.stim.numberreps));
set(S.stim.ampmoddepthbut,'string',num2str(M.stim.ampmoddepth));
set(S.stim.ampmodfreqbut,'string',num2str(M.stim.ampmodfreq));
set(S.stim.samechanbut,'value',M.stim.sameonallchannels);
set(S.stim.waveformbut,'value',M.stim.waveformindex);
set(S.stim.phase1pwbut,'string',num2str(M.stim.phase1pulsewidth));
set(S.stim.phase2pwbut,'string',num2str(M.stim.phase2pulsewidth));
set(S.stim.phase1ampbut,'string',num2str(M.stim.phase1amp));
set(S.stim.phase2ampbut,'string',num2str(M.stim.phase2amp));
set(S.stim.phasegapbut,'string',num2str(M.stim.phasegap));

if isfield(M,'basestim')
      % delay normal stimulus
        S.stim.burstdelay = M.stim.burstdelay;
        % settings for base stimulus
        S.basestim.stim = M.basestim.stim;
        S.basestim.amplitude = M.basestim.amplitude;
        S.basestim.frequency = M.basestim.frequency;
        S.basestim.phase = M.basestim.phase;
        S.basestim.burstdur = M.basestim.burstdur;
        S.basestim.dc = M.basestim.dc;
        S.basestim.ampmoddepth = M.basestim.ampmoddepth;
        S.basestim.ampmodfreq = M.basestim.ampmodfreq;
        S.basestim.ampmodphase = M.basestim.ampmodphase;
        S.basestim.waveformindex = M.basestim.waveformindex;
        S.basestim.phase1pulsewidth = M.basestim.phase1pulsewidth;
        S.basestim.phase2pulsewidth = M.basestim.phase2pulsewidth;
        S.basestim.phase1amp = M.basestim.phase1amp;
        S.basestim.phase2amp = M.basestim.phase2amp;
        S.basestim.phasegap = M.basestim.phasegap;
        S.basestim.waveformindex = M.basestim.waveformindex;
end

set(S.rec.filenamebut,'string',M.macro.localfilename);
set(S.stim.filenamebut,'string',M.macro.localfilename);

NIStim('parseGUI')
NIStim('continuous')
NIStim('waveform')
NIStim('pulsewidth')
NIStim('frequency')
%NIStim('playbackmode')

%--------------------------------------------------------------------------
function NISMmakeInputBox
global S

S.macro.fhdl = findobj('tag','NIStimMacroFig');
if isempty(S.macro.fhdl)
    S.macro.fhdl = figure;
else
    clf(S.macro.fhdl)
end
figure(S.macro.fhdl)
set(S.macro.fhdl,'DeleteFcn','NIStimMacro(''stop'')')
set(S.macro.fhdl,'menubar','none','toolbar','none','numbertitle','off')
set(S.macro.fhdl,'name','NIStimMacro')
set(S.macro.fhdl,'tag','NIStimMacroFig')
set(S.macro.fhdl,'position',[190   552   288   111])
xb = 0.3;
yb = 0.2;

S.macro.stopbut = uicontrol('Parent',S.macro.fhdl,'Style','push','String','Stop','units','normalized','position',[0.1 0.1 xb yb],'callback','NIStimMacro(''stop'')');
S.macro.pausebut = uicontrol('Parent',S.macro.fhdl,'Style','push','String','Pause','units','normalized','position',[0.6 0.1 xb yb],'callback','NIStimMacro(''pause'')');

%--------------------------------------------------------------------------
function NISMstop
global S
S.macro.active = 0;
delete(S.macro.fhdl)

%--------------------------------------------------------------------------
function NISMpause
global S
S.macro.run = 0;
set(S.macro.pausebut,'string','Resume','callback','NIStimMacro(''resume'')')

%--------------------------------------------------------------------------
function NISMresume
global S
S.macro.run = 1;
set(S.macro.pausebut,'string','Pause','callback','NIStimMacro(''pause'')')
