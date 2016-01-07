function P = NIplotMacro(Ain);

accelChInd = [4 5 6];
for n = 1:length(Ain)
    A = NIplotSeqAvg(Ain(n),[2 3],accelChInd,'noplot');
    
    if strcmp(A.seqparametername,'amplitude')
        P.freqVec(n) = A.frequency;
        P.respAmpX(n,:) = A.respAmpX;
        P.respAmpY(n,:) = A.respAmpY;
        P.respAmpZ(n,:) = A.respAmpZ;
        P.respSTDX(n,:) = A.respSTDX;
        P.respSTDY(n,:) = A.respSTDY;
        P.respSTDZ(n,:) = A.respSTDZ;
   
    elseif strcmp(A.seqparametername,'phasegap')
        P.freqVec(n) = A.frequency;
        P.ampVec(n) = A.amplitude;
        P.respAmpX(n,:) = A.respAmpX;
        P.respAmpY(n,:) = A.respAmpY;
        P.respAmpZ(n,:) = A.respAmpZ;
        P.respSTDX(n,:) = A.respSTDX;
        P.respSTDY(n,:) = A.respSTDY;
        P.respSTDZ(n,:) = A.respSTDZ;
   
    elseif strcmp(A.seqparametername,'frequency')
        P.ampVec(n) = A.amplitude;
        P.respAmpX(:,n) = A.respAmpX;
        P.respAmpY(:,n) = A.respAmpY;
        P.respAmpZ(:,n) = A.respAmpZ;
        P.respSTDX(:,n) = A.respSTDX;
        P.respSTDY(:,n) = A.respSTDY;
        P.respSTDZ(:,n) = A.respSTDZ;
    end
end

if strcmp(A.seqparametername,'frequency')
    P.freqVec = A.frequency;
elseif strcmp(A.seqparametername,'amplitude')
    P.ampVec = A.amplitude;
elseif strcmp(A.seqparametername,'phasegap')
    P.ampVec = A.phasegap;
end

% [P.freqVec,sInd] = sort(P.freqVec);
% P.respAmpX = P.respAmpX(sInd,:);
% P.respAmpY = P.respAmpY(sInd,:);
% P.respAmpZ = P.respAmpZ(sInd,:);
% P.respSTDX = P.respSTDX(sInd,:);
% P.respSTDY = P.respSTDY(sInd,:);
% P.respSTDZ = P.respSTDZ(sInd,:);

%------------- Line plots ------------
copo = {'b','r','g','m','k','y','c','b','r','g','m','k','y','c'};
np = size(P.respAmpY,1);
figure('position',[15         127        1342         420]);
subplot(1,3,1)
for n = 1:np
    %plot(P.ampVec,P.respAmpX(n,:),copo{n})
    errorbar(P.ampVec,P.respAmpX(n,:),P.respSTDX(n,:),copo{n})
    hold on
end
xlabel('Amplitude (mA)')
ylabel('Displacement X')
box off
set(gca,'tickdir','out')

subplot(1,3,2)
for n = 1:np
    %plot(P.ampVec,P.respAmpY(n,:),copo{n})
    errorbar(P.ampVec,P.respAmpY(n,:),P.respSTDY(n,:),copo{n})
    hold on
end
xlabel('Amplitude (mA)')
ylabel('Displacement Y')
box off
set(gca,'tickdir','out')


subplot(1,3,3)
for n = 1:np
    %plot(P.ampVec,P.respAmpZ(n,:),copo{n})
    errorbar(P.ampVec,P.respAmpZ(n,:),P.respSTDZ(n,:),copo{n})
    hold on
end
xlabel('Amplitude (mA)')
ylabel('Displacement Z')
box off
set(gca,'tickdir','out')

legend({Ain.macro})




