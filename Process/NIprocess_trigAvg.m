function OUT = NIprocess_trigAvg(command,input)

switch command
    case 'settings'
        OUT = NIPsettings(input);
    case 'plot'
        NIPplot
    case 'data'
        NIPdata(input);
end

%--------------------------------------------------------------------------
function S = NIPsettings(S)

S.proc.rawdata = [];
S.proc.fft = 1;
S.proc.trigavg = 1;
S.proc.currentchind = 1;
S.proc.voltagechind = 2;
S.proc.currentconversion = 10e-3;
S.proc.voltageconversion = 20; %1; % DS5 use 20; SRS 560 Amp use Gain setting
S.proc.triglevel = S.trigger.amp/2;

if S.proc.trigavg == 1
    S.proc.timewindow = [-100 1000];
    S.proc.sampvec = [round(S.proc.timewindow(1)/1e3 * S.ni.rate)+1 : round(S.proc.timewindow(2)/1e3 * S.ni.rate)];
    S.proc.timevec = 1e3*(S.proc.sampvec/S.ni.rate); % [-S.proc.timewindow(1)/1e3+1/S.ni.rate:1/S.ni.rate:S.proc.timewindow(2)/1e3];
    S.proc.sumdata = zeros(length(S.proc.timevec),S.ni.nchin);
    S.proc.avgdata = zeros(length(S.proc.timevec),S.ni.nchin);
    S.proc.navg = 0;
    S.proc.txlims = S.proc.timewindow;
end
if S.proc.fft == 1
    S.proc.nfft = length(S.proc.timevec);
    S.proc.l = length(S.proc.timevec);
    S.proc.fftdata = zeros(S.proc.nfft/2+1,2);
    S.proc.freqvec = S.ni.rate/2*linspace(0,1,S.proc.nfft/2+1);
    S.proc.fxlims = [0 1000];%[min(S.proc.freqvec) max(S.proc.freqvec)];
end

%--------------------------------------------------------------------------
function NIPplot
global S NI

if S.proc.trigavg == 1
    if isfield(S.proc,'a1')
        delete(S.proc.a1)
    end
    S.proc.a1 = axes('parent', S.fig.tab3);
    set(S.proc.a1,'position',[0.1 0.6 0.85 0.35])
    S.proc.p1 = plot(S.proc.timevec,S.proc.avgdata(:,[S.proc.currentchind S.proc.voltagechind]));
    xlabel('Average - Time (ms)')
    ylabel('Amp (V | mA)')
    xlim(S.proc.txlims)
    S.proc.lh = legend([num2str(S.proc.navg) ' reps']);
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
    ylabel('Amp (kOhm)')
    xlim(S.proc.fxlims)
end

%--------------------------------------------------------------------------
function NIPdata(eventData)
global S NI

S.proc.rawdata = [S.proc.rawdata; eventData];
%disp(['Size of raw data = ' num2str(size(S.proc.rawdata))]);

%trigInd = find(S.proc.rawdata(1:end-1,S.ni.nchin)<S.proc.triglevel & S.proc.rawdata(2:end,S.ni.nchin)>S.proc.triglevel);
trigInd = find(S.proc.rawdata(1:end-1,S.trigger.chind)<S.proc.triglevel & S.proc.rawdata(2:end,S.trigger.chind)>S.proc.triglevel);
%disp(['Number of triggers = ' num2str(length(trigInd))]);

lastSampUsed = 0;
for n = 1:length(trigInd)
    %disp(['End trigger sample = ' num2str(S.proc.sampvec(end)+trigInd(n))])
    if trigInd(n)+S.proc.sampvec(end)<length(S.proc.rawdata)
        S.proc.navg = S.proc.navg +1;
        S.proc.sumdata = S.proc.sumdata + S.proc.rawdata(S.proc.sampvec+trigInd(n),:);
        lastSampUsed = trigInd(n);
    end
end
if S.proc.navg>0
    S.proc.avgdata = S.proc.sumdata/S.proc.navg;
end

%disp(['lastSampUsed = ' num2str(lastSampUsed)]);
if isempty(trigInd) % no triggers in this data - throw away
    S.proc.rawdata = S.proc.rawdata(end+S.proc.sampvec(1)-1:end,:);
else
    S.proc.rawdata = S.proc.rawdata(lastSampUsed+1:end,:);
end

S.proc.avgdata(:,S.proc.currentchind) = S.proc.avgdata(:,S.proc.currentchind)*S.proc.currentconversion;
S.proc.avgdata(:,S.proc.voltagechind) = S.proc.avgdata(:,S.proc.voltagechind)*S.proc.voltageconversion;

if S.proc.trigavg == 1
%     for n = [S.proc.currentchind S.proc.voltagechind]%1:S.ni.nchin
%         if n == S.proc.currentchind
%             set(S.proc.p1(n),'ydata',S.proc.avgdata(:,n)*1e3)
%         end
%         set(S.proc.p1(n),'ydata',S.proc.avgdata(:,n))
%     end
    set(S.proc.p1(1),'ydata',S.proc.avgdata(:,S.proc.currentchind)*1e3)
    set(S.proc.p1(2),'ydata',S.proc.avgdata(:,S.proc.voltagechind))
    legend(S.proc.a1,[num2str(S.proc.navg) ' reps'])
    legend boxoff
end
if S.proc.fft == 1
    y = fft(S.proc.avgdata,S.proc.nfft)/S.proc.l;
    S.proc.fftdata = 2*abs(y(1:S.proc.nfft/2+1,:));
    S.proc.fftdata_ang = angle(y(1:S.proc.nfft/2+1,:));
    Z = S.proc.fftdata(:,S.proc.voltagechind)./S.proc.fftdata(:,S.proc.currentchind);
    set(S.proc.p2,'ydata',Z/1e3)
%     set(S.proc.p2(1),'ydata',S.proc.fftdata(:,S.proc.currentchind))
%     set(S.proc.p2(2),'ydata',S.proc.fftdata(:,S.proc.voltagechind))
end