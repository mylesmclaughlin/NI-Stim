list = {'ChirpMC_0p1.proc','ChirpMC_0p2.proc','ChirpMC_0p3.proc','ChirpMC_0p4.proc'};

figure; hold on
copo = {'b','g','r','k'};
for n = 1:length(list)
    P(n) = load(['C:\Users\u0043883\Data\NIData\' list{n}],'-mat');
    plot(P(n).P.freqvec,(P(n).P.fftdata(:,3)./P(n).P.fftdata(:,2)),copo{n})
end
legend(list)
title('Motor  Cortex')

list = {'Chirp_0p07.proc','Chirp_0p09.proc','Chirp_0p12.proc','Chirp_0p2.proc'};

figure; hold on
copo = {'b','g','r','k'};
for n = 1:length(list)
    P(n) = load(['C:\Users\u0043883\Data\NIData\' list{n}],'-mat');
    plot(P(n).P.freqvec,(P(n).P.fftdata(:,3)./P(n).P.fftdata(:,2)),copo{n})
end
legend(list)
title('Whisker')

