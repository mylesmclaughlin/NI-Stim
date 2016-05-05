function OUT = NIprocess_closedLoop(command,input)

switch command
    case 'settings'
        OUT = NIPsettings(input);
    case 'plot'
        NIPplot
    case 'data'
        NIPdata(input)
    case 'startStim'
        NIPstartStim
    case 'stopStim'
        NIPstopStim
end

%--------------------------------------------------------------------------
function S = NIPsettings(S)

S.proc.bufferdur = 20;
S.proc.buffersize = S.ni.rate*S.proc.bufferdur;
S.proc.rawdata = zeros(S.proc.buffersize,S.ni.nchin);
S.proc.chdisp = 2;
S.proc.procdata = zeros(S.proc.buffersize,S.proc.chdisp);   
S.proc.chilabel = {'Current','Displacement'};
S.proc.timeseries = 1;
S.proc.specgram = 1;


if S.proc.timeseries == 1;
    S.proc.dispdur = 2;
    S.proc.dispbuffersize = S.ni.rate*2;
    S.proc.disptimedata = S.proc.procdata(end-S.proc.dispbuffersize+1:end,:); 
    S.proc.timevec = [-S.proc.dispdur:1/S.ni.rate:-1/S.ni.rate];
    S.proc.hxlims = [-S.proc.dispdur 0];
    
    % Kalman params
    S.proc.zerotvec = [0:1/S.ni.rate:S.proc.bufferdur - 1/S.ni.rate];
    S.proc.x0(:,1) = [0; 0];
    S.proc.x0(:,2) = [0; 0];
    S.proc.x0(:,3) = [0; 0];
    
    A = [1 1/S.ni.rate; 0 1];
    B = [0;1/S.ni.rate];
    C = [1 0];
    D = [0];
    ts = 1/S.ni.rate;
    FilterSys = ss(A,B,C,D,ts,'InputName','Acceleration', 'OutputName', 'position');
    Q = [1]; R = [1]; % Variances of the noise: assumption, R berekenen via calibratie
    [S.proc.kalmfdisp, L,P,M] = kalman(FilterSys,Q,R,0);
    
end
if S.proc.specgram == 1;
    S.proc.window = S.ni.rate;
    S.proc.noverlap = S.proc.window*0.5;
    S.proc.nfft = 2^nextpow2(S.proc.window)*2;
    [SP,S.proc.freqvec,S.proc.tvec,P] = spectrogram(S.proc.procdata(:,2), S.proc.window, S.proc.noverlap,S.proc.nfft, S.ni.rate,'yaxis');
    S.proc.psddata = mean(abs(P)');
    [S.proc.termoramp,ind] = max(S.proc.psddata);
    S.proc.termorfreq = S.proc.freqvec(ind);
    S.proc.fxlims = [S.proc.freqvec(1) 50];
    [S.proc.filterb,S.proc.filtera] = butter(2, [3 30]/(S.ni.rate/2), 'bandpass');
end
   
% check for calibration file
if S.accel.donecalib == 1
    S.proc.callibration.g0 = S.accel.callibration.g0;
    S.proc.callibration.g1 = S.accel.callibration.g1;
else
    S.proc.callibration.g0 = [1.65 1.65 1.8];
    S.proc.callibration.g1 = [0.3  0.3  0.28];
end

%--------------------------------------------------------------------------
function NIPstartStim
global S

S.proc.closetheloop = 1;
set(S.proc.closedloopbut,'String','Stop CL','Callback','NIprocess_closedLoop(''stopStim'')','backgroundcolor',[1 0 0])

%--------------------------------------------------------------------------
function NIPstopStim
global S

S.proc.closetheloop = 0;
set(S.proc.closedloopbut,'String','Start CL','Callback','NIprocess_closedLoop(''startStim'')','backgroundcolor',[0 1 0])
if S.stim.stim == 1
    NIStim('stopStim')
end

%--------------------------------------------------------------------------
function NIPplot
global S NI

if S.proc.timeseries == 1
    if isfield(S.proc,'a1')
        delete(S.proc.a1)
    end
    S.proc.a1 = axes('parent', S.fig.tab3);
    set(S.proc.a1,'position',[0.1 0.6 0.85 0.35])
    S.proc.p1 = plot(S.proc.timevec,S.proc.disptimedata);
    xlabel('Time (s)')
    ylabel('Amp (V)')
    xlim(S.proc.hxlims)
    l = legend(S.proc.chilabel);
    legend boxoff
    S.proc.t2 = title('Delay = 0  Phase = 0');
end
if S.proc.specgram == 1
    if isfield(S.proc,'a2')
        delete(S.proc.a2)
    end
    S.proc.a2 = axes('parent', S.fig.tab3);
    set(S.proc.a2,'position',[0.1 0.15 0.85 0.35])
    S.proc.p2 = plot(S.proc.freqvec,S.proc.psddata);
    hold on
    S.proc.p3 = plot(S.proc.termorfreq,S.proc.termoramp,'r.');
    
    xlabel('Frequency (Hz)')
    ylabel('PSD')
    title(['Tremor -  Peak Freq = ' num2str(S.proc.termorfreq) ', Peak Amp = ' num2str(S.proc.termoramp)])
    xlim(S.proc.fxlims)
end

%--------------------------------------------------------------------------
function NIPdata(eventData)
global S NI

S.proc.rawdata = [S.proc.rawdata(S.ni.buffersize+1:end,:); eventData];
if S.proc.timeseries == 1
    
    AccelData = NIPconvert(S.proc.rawdata(:,[4:6]));
    
    % x
    [outx, t, x] = lsim(S.proc.kalmfdisp,AccelData(:,1), S.proc.zerotvec, S.proc.x0(:,1));
    S.proc.x0(:,1) = x(end,:);
    displace(:,1) = x(:,1);
    
    % y
    [outy, t, y] = lsim(S.proc.kalmfdisp,AccelData(:,2), S.proc.zerotvec, S.proc.x0(:,2));
    S.proc.x0(:,2) = y(end,:);
    displace(:,3) = y(:,1);
    
    % z
    [outz, t, z] = lsim(S.proc.kalmfdisp,AccelData(:,3), S.proc.zerotvec, S.proc.x0(:,3));
    S.proc.x0(:,3) = z(end,:);
    displace(:,3) = z(:,1);
    
    S.proc.procdata(:,1) = S.proc.rawdata(:,1);
    S.proc.procdata(:,2) = sqrt(displace(:,1).^2 + displace(:,2).^2 + displace(:,3).^2);
    S.proc.procdata(:,2) = filter(S.proc.filterb,S.proc.filtera,S.proc.procdata(:,2));
    
    S.proc.disptimedata = S.proc.procdata(end-S.proc.dispbuffersize+1:end,:);
    for n = 1:S.proc.chdisp
        set(S.proc.p1(n),'ydata',S.proc.disptimedata(:,n))
    end
    
end

if S.proc.specgram == 1
    [SP,f,t,P] = spectrogram(S.proc.procdata(:,2), S.proc.window, S.proc.noverlap,S.proc.nfft, S.ni.rate,'yaxis');
    S.proc.psddata = mean(abs(P)');
    [S.proc.termoramp,ind] = max(S.proc.psddata);
    S.proc.termorfreq = S.proc.freqvec(ind);
    
    set(S.proc.p2,'ydata',S.proc.psddata)
    set(S.proc.p3,'ydata',S.proc.termoramp)
    set(S.proc.p3,'xdata',S.proc.termorfreq)
    title(['Tremor -  Peak Freq = ' num2str(S.proc.termorfreq) ', Peak Amp = ' num2str(S.proc.termoramp)])
end

if S.proc.closetheloop == 1;
    currentStimFreq = str2num(get(S.stim.freqbut,'string'))
    if S.proc.termorfreq~=currentStimFreq
        set(S.stim.freqbut,'string',num2str(S.proc.termorfreq));
        if S.stim.stim == 0
            NIStim('startStim')
        else
            NIStim('updateStim')
        end
    end
end

%--------------------------------------------------------------------------
function AccelData = NIPconvert(rawdata);
global S

AccelData(:,1) = 9.81*(rawdata(:,1)-S.proc.callibration.g0(1))/S.proc.callibration.g1(1);
AccelData(:,2) = 9.81*(rawdata(:,2)-S.proc.callibration.g0(2))/S.proc.callibration.g1(2);
AccelData(:,3) = 9.81*(rawdata(:,3)-S.proc.callibration.g0(3))/S.proc.callibration.g1(3);

