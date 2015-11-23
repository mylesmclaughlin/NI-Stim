function A = NIimpedanceSeqAvg(A,doplot)

if nargin<2
    doplot = 1;
end
curChInd = 1;
volChInd = 2;
[dum,eInd] = min(abs(A.tvec-A.repDur/1e3));
[dum,sInd] = min(abs(A.tvec-0));
A.tvecShort = A.tvec(sInd:eInd);

nCol = 2;
nRow = ceil(A.nSeq/nCol);
if doplot
    f1 = figure;
    f2 = figure;
    f3 = figure;
end
for n = 1:A.nSeq
    
    
    %A.current(:,n) = 10e-3*(A.allData(n,sInd:eInd,curChInd)-mean(A.allData(n,:,curChInd)));
    %A.voltage(:,n) = A.allData(n,sInd:eInd,volChInd)-mean(A.allData(n,:,volChInd));
    
    A.current(:,n) = 10e-3*(A.avgData(n,sInd:eInd,curChInd)-mean(A.avgData(n,:,curChInd)));
    A.voltage(:,n) = A.avgData(n,sInd:eInd,volChInd)-mean(A.avgData(n,:,volChInd));
    
    
    A.impedance(:,n) = A.voltage(:,n)./A.current(:,n);
    [A.xcorr(:,n),A.lags] = xcorr(A.current(:,n),A.voltage(:,n));
    A.impedanceAmpTD(n) = mean(A.impedance(:,n));
    %A.impedanceAmp(n)  = max(mfft(A.voltage(:,n),A.fs))/max(mfft(A.current(:,n),A.fs));
    if strcmp(A.seqparametername,'frequency')
        %A.impedanceAmp(n)  = max(mfft(A.voltage(:,n),A.fs))/(A.amplitude*1e-3);
        [Vamp,VampInd,f,Ymag] = fftM(A.voltage(:,n),A.fs,A.frequency(n));
        A.impedanceAmp(n)  = Vamp/(A.amplitude*1e-3);
    elseif strcmp(A.seqparametername,'amplitude')
        A.impedanceAmp(n)  = fftM(A.voltage(:,n),A.fs,A.frequency)/(A.amplitude(n)*1e-3);
    elseif strcmp(A.seqparametername,'ampmoddepth')
        A.impedanceAmp(n)  = fftM(A.voltage(:,n),A.fs,A.frequency)/(A.amplitude*1e-3);
    end
    [dum,A.lagInd(n)] = max(A.xcorr(:,n));
    A.impedanceDelay(n) = A.lags(A.lagInd(n))/A.fs;
    
    if doplot
        figure(f1)
        subplot(nRow,nCol,n)
        ydata = smooth(A.impedance(:,n),500);
        plot(A.tvecShort*1e3,ydata)
        xlim([A.tvecShort(1) A.tvecShort(end)]*1e3)
        
        ymax(n) = max(ydata(1000:end-1000));
        ymin(n) = min(ydata(1000:end-1000));
        
        if n==A.nSeq
            xlabel('Time (ms)')
        end
        
        figure(f2)
        subplot(nRow,nCol,n)
        plot(A.lags,A.xcorr(:,n))
        hold on
        plot(A.lags(A.lagInd(n)),A.xcorr(A.lagInd(n),n),'ro')
        
        figure(f3)
        subplot(nRow,nCol,n)
        plot(f,Ymag)
        hold on
        plot(f(VampInd),Ymag(VampInd),'r.')
    end
end

A.voltageRMS = rms(A.voltage);
A.currentRMS = rms(A.current);

%yl = [0 max(ymax)];

% for n = 1:A.nSeq
%     figure(f1)
%     subplot(nRow,nCol,n)
%     ylim(yl)
%
% end

%-------------------------------------------------------------
function [Amp,AmpInd,f,Ymag] = fftM(y,fs,Fc)

% append zeros to make length a multiple of Fc
Pc = 1/Fc;
yT = (length(y))/fs;
Ratio = yT/Pc;
ZerosT = ceil(Ratio)-Ratio;
ZerosSamp = round(ZerosT*Pc*fs);
NFFT = length(y)+ZerosSamp;
HNFFT = round(NFFT/2);
y = [y; zeros(ZerosSamp,1)];

%NFFT = length(y);
Y = fft(y,NFFT)/NFFT;
Ymag = 2*abs(Y(1:HNFFT+1));
f = fs/2*linspace(0,1,HNFFT+1);

[dum,AmpInd] = min(abs(Fc-f));
Amp = Ymag(AmpInd);

