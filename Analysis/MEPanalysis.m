
function [A,P] = MEPanalysis(fileName, filePath, tacsFreq, TrainPulseDuration)


A = NImakeSeqAvg(fileName, filePath, 4:6);

ch =4;
plotShift = [1 0.05 0.05];
RepNumber = length(squeeze(A(1).allData(:,1,1,1)));
colorIn = ['b','r','g','m','k'];


TrainPulseDurationSamples = round(TrainPulseDuration/1000 * A(1).fs + 0.2*A(1).fs); % 100 ms before and after the train
PulseNumber = round((TrainPulseDuration/1000)*A(1).frequency);


%Plotting
for j =1:length(A(1).amplitude)
    for i = 1:length(A)
        figure
        shift = round(A(i).phase * 1/tacsFreq *A(i).fs+1);
        plot(squeeze(A(i).allData(1,j,shift:TrainPulseDurationSamples+shift,ch)))
        hold on
        for  k =2:RepNumber
            plot(squeeze(A(i).allData(k,j,shift:TrainPulseDurationSamples+shift,ch)- (k-1)*plotShift(ch-3)),colorIn(k))
        end
        for z = 0:PulseNumber-1
            PulseIndex = 0.1*A(1).fs+1/(A(1).frequency)*z*A(i).fs;
            minData = min(A(i).allData(RepNumber,j,shift:TrainPulseDurationSamples+shift,ch)-(RepNumber-1)*plotShift(ch-3));
            maxData = max(A(i).allData(1,j,TrainPulseDurationSamples+shift,ch));
            line([PulseIndex PulseIndex],[minData maxData],'LineStyle','--');
        end
        axis ('tight')
    end
end



% Area calculation 50 ms after last Pulse

wo = 50/(A(1).fs/2);  bw = wo/10;
[b1,a1] = iirnotch(wo,bw); % 50 hz


%[b,a] = butter(2,10/(A(1).fs/2),'low');
[b,a] = butter(2,[30 1000]/(A(1).fs/2),'bandpass');
%Window = [-50 200]/1000; % 
Window = [0.002 0.05];

for j = 1:length(A);
    LastPulseIndex(j) = round(0.1*A(j).fs + (PulseNumber-1) * 1/A(j).frequency * A(j).fs + A(j).phase * 1/tacsFreq *A(i).fs); %1/(A(1).frequency)*j*A(i).fs;
    
    for i = 1:length(A(1).amplitude)
        for z = 1:RepNumber
            data = filter(b1,a1,A(j).allData(z,i,:,ch));
            data = filter(b,a,data);
            data1 = data(LastPulseIndex(j)+Window(1)*A(j).fs:LastPulseIndex(j)+Window(2)*A(j).fs);
            A(j).AreaMEP(i,z) = sum(abs(data1));
            [pks,locs] = findpeaks(squeeze(data1),'MinPeakHeight',0.05);
%             A(j).pks(i,z,:) = pks(1:3);
%             A(j).locs(i,z,:) = locs(1:3);
%             A(j).locs(i,z,:) = A(j).locs(i,z,:) + LastPulseIndex(j)+Window(1)*A(j).fs;
        end
    end
%     figure
%     plot(squeeze(A(j).allData(RepNumber,length(A(1).amplitude),:,4)))
%     hold on
%     plot(squeeze(data),'g')
%     line([LastPulseIndex(j)+Window(1)*A(j).fs LastPulseIndex(j)+Window(1)*A(j).fs],[min(A(j).allData(2,2,:,4)) max(A(j).allData(2,2,:,4))],'color','r','LineStyle','--');
%     line([LastPulseIndex(j)+Window(2)*A(j).fs LastPulseIndex(j)+Window(2)*A(j).fs],[min(A(j).allData(2,2,:,4)) max(A(j).allData(2,2,:,4))],'color','r','LineStyle','--');
    A(j).MeanMEP = mean(A(j).AreaMEP.');
    A(j).StdMEP = std(A(j).AreaMEP',0,1);
end

for j = 1:length(A)
    for i = 1:length(A(1).amplitude)
        a = A(j).AreaMEP(i,:);
        a(a > mean(a)+2*std(a) | a < mean(a)-2*std(a))=[];
        A(j).AreaMeanNoOutliers(:,i) = mean(a.');
        A(j).AreaSTDNoOutliers(:,i) = std(a.',0,1);
    end
end

% plot

for j = 1:length(A)
    MeanMEPvec(j,:) = A(j).MeanMEP;
    StdMEPvec(j,:) = A(j).StdMEP;
end



% figure
% 
% for i=1:length(MeanMEPvec(1,:))
%     if i==1
%         errorbar(MeanMEPvec(:,1),StdMEPvec(:,i))
%         hold on
%     else
%         errorbar(MeanMEPvec(:,i),StdMEPvec(:,i),colorIn(i))
%     end
% end

%%
phaseplot = 1;
sortFreq = 1;
Ain =A;
if sortFreq
    for n = 1:length(Ain)
        freq(n) = Ain(n).frequency;
    end
    [s,sind] = sort(freq);
    Ain = Ain(sind);
end

accelChInd = [4 5 6];
for n = 1:length(Ain)
    A1 = NIplotSeqAvg(Ain(n),[2 3],accelChInd,'noplot');
    %A = NIplotSeqAvg(Ain(n),[2 3],accelChInd,'allPlot');
    if strcmp(A1.seqparametername,'amplitude')
        P.freqVec(n) = A1.frequency;
        [P.respAmpX(n,:),ind] =  max(A1.avgData(:,:,4)');
        for i = 1:length(ind)
            P.respSTDX(n,i) = A1.stdData(i,ind(i),4);
        end
        
        %P.respAmpX(n,:) = A.respAmpX;
        P.respAmpY(n,:) = A1.respAmpY;
        P.respAmpZ(n,:) = A1.respAmpZ;
        %P.respSTDX(n,:) = A.respSTDX;
        P.respSTDY(n,:) = A1.respSTDY;
        P.respSTDZ(n,:) = A1.respSTDZ;
        P.respAllX(n,:,:) = A1.respAllX;
        P.respAllY(n,:,:) = A1.respAllY;
        P.respAllZ(n,:,:) = A1.respAllZ;
        P.StdMEP (n,:) = A1.StdMEP;
        P.MeanMEP (n,:) = A1.MeanMEP;
        P.MeanMEPvec(n,:) = A1.MeanMEP;
        P.StdMEPvec(n,:) = A1.StdMEP;
    elseif strcmp(A1.seqparametername,'phasegap')
        P.freqVec(n) = A1.frequency;
        P.ampVec(n) = A1.amplitude;
        P.respAmpX(n,:) = A1.respAmpX;
        P.respAmpY(n,:) = A1.respAmpY;
        P.respAmpZ(n,:) = A1.respAmpZ;
        P.respSTDX(n,:) = A1.respSTDX;
        P.respSTDY(n,:) = A1.respSTDY;
        P.respSTDZ(n,:) = A1.respSTDZ;
        P.respAllX(n,:,:) = A1.respAllX;
        P.respAllY(n,:,:) = A1.respAllY;
        P.respAllZ(n,:,:) = A1.respAllZ;
        P.StdMEP (n,:) = A1.StdMEP;
        P.MeanMEP (n,:) = A1.MeanMEP;
        P.MeanMEPvec(n,:) = A1.MeanMEP;
        P.StdMEPvec(n,:) = A1.StdMEP;
    elseif strcmp(A1.seqparametername,'frequency')
        P.ampVec(n) = A1.amplitude;
        P.respAmpX(:,n) = A1.respAmpX;
        P.respAmpY(:,n) = A1.respAmpY;
        P.respAmpZ(:,n) = A1.respAmpZ;
        P.respSTDX(:,n) = A1.respSTDX;
        P.respSTDY(:,n) = A1.respSTDY;
        P.respSTDZ(:,n) = A1.respSTDZ;
        P.respAllX(n,:,:) = A1.respAllX;
        P.respAllY(n,:,:) = A1.respAllY;
        P.respAllZ(n,:,:) = A1.respAllZ;
        P.StdMEP (n,:) = A1.StdMEP;
        P.MeanMEP (n,:) = A1.MeanMEP;
        P.MeanMEPvec(n,:) = A1.MeanMEP;
        P.StdMEPvec(n,:) = A1.StdMEP;
    end
end

selectAxis(1) = max(P.respAmpX(:));
selectAxis(2) = max(P.respAmpY(:));
selectAxis(3) = max(P.respAmpZ(:));
[dum,ind] = max(selectAxis);
if ind == 1
    P.respAmp = P.respAmpX;
    P.respAll = P.respAllX;
    P.respSTD = P.respSTDX;
elseif ind == 2
    P.respAmp = P.respAmpY;
    P.respAll = P.respAllY;
    P.respSTD = P.respSTDY;
elseif ind == 3
    P.respAmp = P.respAmpZ;
    P.respAll = P.respAllZ;
    P.respSTD = P.respSTDZ;
end


if strcmp(A1.seqparametername,'frequency')
    P.freqVec = A1.frequency;
elseif strcmp(A1.seqparametername,'amplitude')
    P.ampVec = A1.amplitude;
elseif strcmp(A1.seqparametername,'phasegap')
    P.ampVec = A1.phasegap;
end



%------------- Line plots ------------
copo = {'b','r','g','m','k','y','c','b','r','g','m','k','y','c'};
np = size(P.respAmpY,1);
nq = size(P.respAmpY,2);

if phaseplot
    
    P.phase = [Ain.phase];
    P.subthresholdfrequency = 2;
    copo = {'b','r','g','m','k','y'};
    figure('position',[15         127        1342         420])
    np = length(P.respAmpX(1,:));
    subplot(1,3,1)
    for i = 1:np
        
        %errorbar(P.phase,P.respAmpX(:,i),P.respSTDX(:,i),[copo{i} '.-'])
        errorbar(P.phase,P.MeanMEPvec(:,i),P.StdMEPvec(:,i),[copo{i} '.-'])
        
        hold on
    end
    xlabel('Phase')
    ylabel('AUC (a.a)')
    title('MEP')
    
    subplot(1,3,2)
    for i = 1:np
        errorbar(P.phase,P.respAmpY(:,i),P.respSTDY(:,i),[copo{i} '.-'])
        hold on
    end
    xlabel('Phase')
    title(['Y-axis' num2str(P.subthresholdfrequency) ' Hz subthreshold stim'])
    
    subplot(1,3,3)
    for i = 1:np
        errorbar(P.phase,P.respAmpZ(:,i),P.respSTDZ(:,i),[copo{i} '.-'])
        hold on
    end
    xlabel('Phase')
    title(['Z-axis' num2str(P.subthresholdfrequency) ' Hz subthreshold stim'])
    
end

