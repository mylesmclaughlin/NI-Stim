function [A]=plotAvgMEP(A)

Ampl = 3;% to plot

ch = 4;
AcclCh = [5];%[4 5 6];
colorInd = {'b','g','r','m','k','b','g','r','m','k','b','g','r','m','k','b','g','r','m','k'};
sampWindow = [A(1).fs*0.1 A(1).fs*0.3];
sampWindowMEP = [A(1).fs*0.115 A(1).fs*0.15];
[b,a] = butter(2,[30 1000]/(A(1).fs/2),'bandpass');

% for j = 1:length(A(1).amplitude)
%     for i =1:length(A)
%         avgData = A(i).avgData(j,:,ch);
%         for k =1:length(A(1).allData(:,1,1,ch))+1
%             if k<(length(A(1).allData(:,1,1,ch))+1)
%                 data = filter(b,a,squeeze(A(i).allData(k,j,:,ch)));
%             end
%             if (i==1 && k==1)
%                 figure
%                 plot(data)
%                 hold on
%             else
%                 if (k==length(A(1).allData(:,1,1,ch))+1)
%                     %plot(avgData, 'color', [0.5 0.5 0.5],'LineWidth',4,'LineStyle','--');
%
%                 else
%                     plot(data,'color',colorInd{k})
%                 end
%             end
%             %[A(i).MaxAvgMEP(j),ind] = max(data);
%             %A(i).MaxMEPSTD(j) = std(squeeze(A(i).allData(:,j,ind,ch)));
%         end
%         legend('Rep1','Rep2','Rep3','Rep4','Rep5')%,'Average');
%     end
% end


for j = Ampl:Ampl %1:length(A(1).amplitude)
    for i =1:length(A)
        figure;
        subplot(1,2,1)
        shift = round(A(i).phase/A(i).baseStimFreq*A(i).fs);
        wind = [sampWindowMEP(1)+shift sampWindowMEP(2)+shift];
        avgData = A(i).avgData(j,:,ch)-mean(A(i).avgData(j,:,ch));
        avgData = avgData(wind(1):wind(2));
        for k =1:length(A(1).allData(:,1,1,ch))+1
            if k<(length(A(1).allData(:,1,1,ch))+1)
                data = filter(b,a,squeeze(A(i).allData(k,j,:,ch)));
                data = data (wind(1):wind(2));
            end
            if k==1
                plot(data)
                hold on
            else
                if (k==length(A(1).allData(:,1,1,ch))+1)
                    %plot(avgData, 'color', [0.5 0.5 0.5],'LineWidth',4,'LineStyle','--');
                    
                else
                    plot(data,'color',colorInd{k})
                end
            end
            %[A(i).MaxAvgMEP(j),ind] = max(data);
            %A(i).MaxMEPSTD(j) = std(squeeze(A(i).allData(:,j,ind,ch)));
        end
        legend('Rep1','Rep2','Rep3','Rep4','Rep5')%,'Average');
        axis('tight')
        
        
        subplot(1,2,2)
        wind = [sampWindow(1)+shift sampWindow(1)+6000+shift];
        avgData = A(i).avgData(j,:,AcclCh)-mean(A(i).avgData(j,:,AcclCh));
        avgData = avgData(wind(1):wind(2));
        for k =1:length(A(1).allData(:,1,1,AcclCh))+1
            if k<(length(A(1).allData(:,1,1,AcclCh))+1)
                data = filter(b,a,squeeze(A(i).allData(k,j,:,AcclCh)));
                data = data (wind(1):wind(2));
            end
            if k==1
                plot(data)
                hold on
            else
                if (k==length(A(1).allData(:,1,1,AcclCh))+1)
                    %plot(avgData, 'color', [0.5 0.5 0.5],'LineWidth',2,'LineStyle','--');
                    
                    
                else
                    plot(data,'color',colorInd{k})
                end
            end
            %[A(i).MaxAvgMEP(j),ind] = max(data);
            %A(i).MaxMEPSTD(j) = std(squeeze(A(i).allData(:,j,ind,ch)));
        end
        legend('Rep1','Rep2','Rep3','Rep4','Rep5')%,'Average');
        axis('tight');
    end
end


%[b,a] = butter(2,[1 500]/(A(1).fs/2),'bandpass');





% phases = [A(:).phase];
% MaxAvgMEP = [A(:).MaxAvgMEP];
% MaxMEPSTD = [A(:).MaxMEPSTD];
% L = length(MaxAvgMEP)/length(A(1).amplitude);
% MaxAvgMEP = reshape(MaxAvgMEP,[length(A(1).amplitude),L])';
% MaxMEPSTD = reshape(MaxMEPSTD,[length(A(1).amplitude),L])';
%
% figure('position',[15         127        1342         420]);
% for j = 1:length(A(1).amplitude)
%     errorbar(phases,MaxAvgMEP(:,j),MaxMEPSTD(:,j),colorInd{j})
%     hold on
% end
%
% xlabel('Phase')
% ylabel('MEP max')
% box off
% set(gca,'tickdir','out')