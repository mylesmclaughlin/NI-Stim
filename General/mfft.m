function [Ymag,f] = mfft(y,fs,xrange,doplot)

if nargin<3
    xrange = [0 fs/2];
end
if nargin<4
    doplot = 1;
end

L = length(y);
%NFFT = 2^nextpow2(L); % Next power of 2 from length of y
NFFT = L;
Y = fft(y,NFFT)/L;
Ymag = 2*abs(Y(1:round(NFFT/2)+1));
YmagdB = mag2db(Ymag);
f = fs/2*linspace(0,1,round(NFFT/2)+1);

if doplot
    
    [dum ind(1)] = min(abs(f-xrange(1)));
    [dum ind(2)] = min(abs(f-xrange(2)));
    YmagdBnorm = YmagdB - max(YmagdB(ind(1):ind(2)));

    % Plot single-sided amplitude spectrum.
    fh = findobj('name','fft plot');
    copo = {'b','g','r','c','m','k','y','b','g','r','c','m','k','y'};
    if isempty(fh)
        fh = figure('Name','fft plot','tag','1');
        n = 1;
    else
        figure(fh)
        n = str2num(get(fh,'tag'))+1;
        set(fh,'tag',num2str(n));
    end
    
    plot(f,Ymag,copo{n})
    hold on
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    %ylabel('|Y(f)| dB')
    %set(gca,'xscale','log')
    xlim(xrange)
end