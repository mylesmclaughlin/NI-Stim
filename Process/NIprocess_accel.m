function OUT = NIprocess_accel(command,input)

switch command
    case 'settings'
        OUT = NIPsettings(input);
    case 'plot'
        NIPplot
    case 'data'
        NIPdata(input)
end

%--------------------------------------------------------------------------
function S = NIPsettings(S)

S.proc.dur = 2;
S.proc.buffersize = S.ni.rate*S.proc.dur;
S.proc.rawdata = zeros(S.proc.buffersize,S.ni.nchin);
S.proc.chdisp = 5;
S.proc.chilabel = {'Current','X','Y','Z','Current - X'}
S.proc.fft = 1;
S.proc.hilbert = 1;
if S.proc.hilbert == 1
    S.proc.hilbertdata = zeros(S.proc.buffersize,S.proc.chdisp);
    S.proc.timevec = [-S.proc.dur:1/S.ni.rate:-1/S.ni.rate];
    S.proc.hxlims = [-S.proc.dur 0];
    
    % Kalman params
    S.proc.zerotvec = [0:1/S.ni.rate:S.proc.dur - 1/S.ni.rate];
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
if S.proc.fft == 1
    S.proc.nfft = S.proc.buffersize;
    S.proc.l = S.proc.buffersize;
    S.proc.fftdata = zeros(S.proc.nfft/2+1,S.proc.chdisp-1);
    S.proc.navg = 10;
    S.proc.fftdatabuffer = zeros(S.proc.navg,S.proc.nfft/2+1,S.proc.chdisp-1);
    S.proc.fftdatamean = zeros(size(S.proc.fftdata));
    S.proc.fftdatabufferpos = 0;
    S.proc.freqvec = S.ni.rate/2*linspace(0,1,S.proc.nfft/2+1);
    S.proc.fxlims = [1 25];%[min(S.proc.freqvec) max(S.proc.freqvec)];
end

% check for calibration file
S.proc.callibration.g0 = [1.65 1.65 1.8];
S.proc.callibration.g1 = [0.3  0.3  0.28];

%--------------------------------------------------------------------------
function NIPplot
global S NI

if S.proc.hilbert == 1
    if isfield(S.proc,'a1')
        delete(S.proc.a1)
    end
    S.proc.a1 = axes('parent', S.fig.tab3);
    set(S.proc.a1,'position',[0.1 0.6 0.85 0.35])
    S.proc.p1 = plot(S.proc.timevec,S.proc.hilbertdata);
    xlabel('Hilbert - Time (s)')
    ylabel('Phase (cycles)')
    xlim(S.proc.hxlims)
    l = legend(S.proc.chilabel);
    legend boxoff
end
if S.proc.fft == 1
    if isfield(S.proc,'a2')
        delete(S.proc.a2)
    end
    S.proc.a2 = axes('parent', S.fig.tab3);
    set(S.proc.a2,'position',[0.1 0.15 0.85 0.35])
    S.proc.p2 = plot(S.proc.freqvec,S.proc.fftdata);
    xlabel('FFT - Frequency (Hz)')
    ylabel('Amp (V)')
    xlim(S.proc.fxlims)
end
% ylim(S.rec.ylims)
% ylabel('Voltage (V)')

%--------------------------------------------------------------------------
function NIPdata(eventData)
global S NI

S.proc.rawdata = [S.proc.rawdata(S.ni.buffersize+1:end,:); eventData];
if S.proc.hilbert == 1
    %fdata = filter(S.accel.filterb,S.accel.filtera,S.proc.rawdata);
    fdata = S.proc.rawdata(:,[1 4:6]);
    
    AccelData = NIPconvert(fdata(:,2:4));
    
    % x
    [outx, t, x] = lsim(S.proc.kalmfdisp,AccelData(:,1), S.proc.zerotvec, S.proc.x0(:,1));
    S.proc.x0(:,1) = x(end,:);
    fdata(:,2) = x(:,1);
    
    % y
    [outy, t, y] = lsim(S.proc.kalmfdisp,AccelData(:,2), S.proc.zerotvec, S.proc.x0(:,2));
    S.proc.x0(:,2) = y(end,:);
    fdata(:,3) = y(:,1);
    
    % z
    [outz, t, z] = lsim(S.proc.kalmfdisp,AccelData(:,3), S.proc.zerotvec, S.proc.x0(:,3));
    S.proc.x0(:,3) = z(end,:);
    fdata(:,4) = z(:,1);
    
    
    %fdata = S.proc.rawdata;
    y = hilbert(fdata);
    S.proc.hilbertdata(:,1:S.proc.chdisp-1) = angle(y);
    S.proc.hilbertdata(:,S.proc.chdisp) = (S.proc.hilbertdata(:,1) - S.proc.hilbertdata(:,2));
    S.proc.hilbertdata = S.proc.hilbertdata/(2*pi);
    %S.proc.hilbertdata(:,1) = abs(y);
    for n = 1:S.proc.chdisp
        set(S.proc.p1(n),'ydata',S.proc.hilbertdata(:,n))
    end
end
if S.proc.fft == 1
    
    %y = fft(S.proc.rawdata,S.proc.nfft)/S.proc.l;
    y = fft(fdata,S.proc.nfft)/S.proc.l;
    S.proc.fftdata = 2*abs(y(1:S.proc.nfft/2+1,:));
    
    % store some ffts for averaging
    S.proc.fftdatabufferpos = S.proc.fftdatabufferpos +1;
    if S.proc.fftdatabufferpos > S.proc.navg
        S.proc.fftdatabufferpos = 1;
    end
    S.proc.fftdatabuffer(S.proc.fftdatabufferpos,:,:) = S.proc.fftdata;
    S.proc.fftdatamean = squeeze(mean(S.proc.fftdatabuffer,1));
    for n = 1:S.proc.chdisp-1
        set(S.proc.p2(n),'ydata',S.proc.fftdata(:,n))
    end
end

%--------------------------------------------------------------------------
function AccelData = NIPconvert(rawdata);
global S

AccelData(:,1) = 9.81*(rawdata(:,1)-S.proc.callibration.g0(1))/S.proc.callibration.g1(1);
AccelData(:,2) = 9.81*(rawdata(:,2)-S.proc.callibration.g0(2))/S.proc.callibration.g1(2);
AccelData(:,3) = 9.81*(rawdata(:,3)-S.proc.callibration.g0(3))/S.proc.callibration.g1(3);
