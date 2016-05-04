function NIsinglePhasePlot(A,freq)

if nargin<2
    freq = 2;
end

phase = freq *A.seqparametervalues/1000;
accelInd = [4 5 6];
figure('position',[159 295 1304 420])
for n = 1:3
    val = [];
    ind = [];
    sd = [];
    [val,ind] = max(abs(A.avgDisp(:,:,accelInd(n)))');
    for m = 1:length(ind)
        sd(m) = A.stdDisp(m,ind(m),accelInd(n));
    end
    subplot(1,3,n)
    errorbar(phase,val,sd)
    ylabel('Peak Limb Displacement')
    xlabel('Phase')
    xlim([min(phase) max(phase)])
end
    