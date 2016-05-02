function [D,S] = quickSpecComp(fileName)

load(fileName)
D.tvec = [1:length(D.data)]/D.fs;

[b,a] = butter(2,[2 500]/(D.fs/2),'bandpass');
%[b,a] = butter(2,500/(D.fs/2),'low');
D.data(:,1) = filter(b,a,D.data(:,1));
D.data(:,2) = filter(b,a,D.data(:,2));

D.stimOn = D.stimOn(1:length(D.data));
indOn = find(D.stimOn==1);
indOff = find(D.stimOn==0);
sigOn = D.data(indOn,:);
sigOff = D.data(indOff,:);

window = D.fs;
noverlap = window*0.75;
nfft = 2^nextpow2(window)*2;
flim = [0 100];


figure

subplot(2,2,1)
[S.BST,f,t,P.BST] = spectrogram(D.data(:,1),window,noverlap,nfft,D.fs);
s=surf(t,f, 10*log10(abs(P.BST)), 'linestyle','none');
view([0 90])
hold on
plot(D.tvec,flim(end)/2*D.stimOn,'k');
ylim(flim)
xlim([min(t) max(t)])
colormap hsv
xlabel('Time (S)')
ylabel('Frequency (Hz)')
title('BST: Spectrogram')

subplot(2,2,2)
[S.BSTOff,f,t,P.BSTOff] = spectrogram(sigOff(:,1),window,noverlap,nfft,D.fs);
if ~isempty(sigOn)
    [S.BSTOn,f,t,P.BSTOn] = spectrogram(sigOn(:,1),window,noverlap,nfft,D.fs);
else
    S.BSTOn = [];
    P.BSTOn = [];
end

plot(f,mean(abs(S.BSTOff)'),'b')
hold on
plot(f,mean(abs(S.BSTOn)'),'r')
hold on
xlim(flim)
ylabel('Magnitude')
xlabel('Frequency (Hz)')
title('BST: PSD')
legend('BST: Stim Off','BST: Stim On')

subplot(2,2,3)
[S.NonBST,f,t,P.NonBST] = spectrogram(D.data(:,2),window,noverlap,nfft,D.fs);
s=surf(t,f, 10*log10(abs(P.NonBST)), 'linestyle','none');
view([0 90])
hold on
plot(D.tvec,flim(end)/2*D.stimOn,'k');
ylim(flim)
xlim([min(t) max(t)])
colormap hsv
xlabel('Time (S)')
ylabel('Frequency (Hz)')
title('Non-BST: Spectrogram')

subplot(2,2,4)
[S.NonBSTOff,f,t,P.NonBSTOff] = spectrogram(sigOff(:,2),window,noverlap,nfft,D.fs);
if ~isempty(sigOn(:,2))
    [S.NonBSTOn,f,t,P.NonBSTOn] = spectrogram(sigOn(:,2),window,noverlap,nfft,D.fs);
else
    S.NonBSTOn = [];
    P.NonBSTOn = [];
end
plot(f,mean(abs(S.NonBSTOff)'),'b')
hold on
plot(f,mean(abs(S.NonBSTOn)'),'r')
hold on
xlim(flim)
ylabel('Magnitude')
xlabel('Frequency (Hz)')
title('Non-BSD: PSD')
legend('Non-BST: Stim Off','Non-BST: Stim On')



