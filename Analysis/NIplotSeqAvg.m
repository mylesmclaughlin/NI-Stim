function A = NIplotSeqAvg(A,chInd,accelInd,plotString)

if nargin<2
    chInd = A.nCh;
end
if nargin<3
    accelInd = [];
end
if nargin<4
    plotString = 'allPlot';
end


nCol = 2;
nRow = ceil(A.nSeq/nCol);

for i = 1:length(chInd)
    if strcmpi(plotString,'allPlot')
        figure('name',['Ch ' num2str(chInd(i))],'tag',['NIplotSeqAvg1_' num2str(i)])
    end
    ymax = [];
    ymin = [];
    for n = 1:A.nSeq
        ydata = squeeze(A.avgData(n,:,chInd(i)));
        if strcmpi(plotString,'allPlot')
            subplot(nRow,nCol,n)
            plot(A.tvec*1e3,ydata)
            xlim([A.tvec(1) A.tvec(end)]*1e3)
            if n==A.nSeq
                xlabel('Time (ms)')
            end
        end
        ymax(n) = max(ydata);
        ymin(n) = min(ydata);
        
    end
    yl = [min(ymin) max(ymax)];
    
    if strcmpi(plotString,'allPlot')
        for n = 1:A.nSeq
            subplot(nRow,nCol,n)
            ylim(yl)
        end
    end
end

if ~isempty(accelInd)
    accelString = {'X','Y','Z'};
    bandpass = [1 500];
    filterorder = 2;
    [filterb filtera] = butter(filterorder,bandpass/(A.fs/2),'bandpass');
    
    for n = 1:A.nSeq
        dataX = squeeze(A.avgData(n,:,accelInd(1)));
        dataY = squeeze(A.avgData(n,:,accelInd(2)));
        dataZ = squeeze(A.avgData(n,:,accelInd(3)));
        data = sqrt(dataX.^ + dataY.^2 + dataZ.^2);
        zInd = find(A.tvec<0);
        stimStartInd = zInd(end)+1;
        calib = mean(data(zInd));
        data = data - calib;
        respAccel = filter(filterb,filtera,data);
        data = cumtrapz(respAccel);
        respVel = filter(filterb,filtera,data);
        respDisp = cumtrapz(respVel);
        
        A.respAccel(n,:) = respAccel;
        A.respVel(n,:) = respVel;
        A.respDisp(n,:) = respDisp;
        
    end
    
    for i = 1:length(accelInd)
        if strcmpi(plotString,'allPlot')
            figure('name',['Acceleration: ' accelString{i} ' axis'],'tag',['NIplotSeqAvg2_' num2str(i)])
        end
        ymax = [];
        ymin = [];
        for n = 1:A.nSeq
            
            
            zInd = find(A.tvec<0);
            stimStartInd = zInd(end)+1;
            
            eval(['allData = squeeze(A.allData' num2str(n) '(:,:,accelInd(i)));'])
            
            for j = 1:length(allData(:,1))
                calib = mean(allData(j,zInd));
                allData(j,:) = allData(j,:) - calib;
                %allData(j,artInd) = 0;
                allData(j,:) = filter(filterb,filtera,allData(j,:));
                allData(j,:) = cumtrapz(allData(j,:)); % velocity
                allData(j,:) = filter(filterb,filtera,allData(j,:));
                allData(j,:) = cumtrapz(allData(j,:)); % displacement
                sumData(j) = sum(abs(allData(j,stimStartInd:end)));
            end
            
            ydata = mean(allData);
            ydata_std = std(allData);
            ydataE = allData(end,:);
            ydataS = allData(1,:);
            
            ymax(n) = max(ydata);
            ymin(n) = min(ydata);
            [respAmp(n),ind] = max(abs(ydata(stimStartInd:end)));
            respSTD(n) = ydata_std(ind);
            respAll(n,:) = max(abs(allData(:,stimStartInd:end))');
           
%             respAmp(n) = mean(sumData);
%             respSTD(n) = std(sumData);
            respAmpE(n) = max(abs(ydataE(stimStartInd:end)));
            respAmpS(n) = max(abs(ydataS(stimStartInd:end)));
            
            if strcmpi(plotString,'allPlot')
                subplot(nRow,nCol,n)
                plot(A.tvec*1e3,ydata)
                hold on
                plot(A.tvec*1e3,ydataE,'r')
                plot(A.tvec*1e3,ydataS,'g')
                xlim([A.tvec(1) A.tvec(end)]*1e3)
                if n==A.nSeq
                    xlabel('Time (ms)')
                end
            end
            
            
        end
        yl = [min(ymin) max(ymax)];
        
        if strcmpi(plotString,'allPlot')
            for n = 1:A.nSeq
                subplot(nRow,nCol,n)
                ylim(yl)
            end
        end
        
        if strcmpi(plotString,'allPlot')
            figure('name',['Input-Output Function - Acceleration: ' accelString{i} ' axis'],'numbertitle','off','tag',['NIplotSeqAvg3_' num2str(i)])
            plot(A.seqparametervalues,respAmp,'.-')
            hold on
            plot(A.seqparametervalues,respAmpE,'r.-')
            plot(A.seqparametervalues,respAmpS,'g.-')
            xlabel(['INPUT:' A.seqparametername])
            ylabel(['OUTPUT: Response Amplitude'])
            title(['Input-Output Function - Acceleration: ' accelString{i} ' axis'])
        end
        eval(['A.respAmp' accelString{i} ' = respAmp;']);
        eval(['A.respAll' accelString{i} ' = respAll;']);
        eval(['A.respSTD' accelString{i} ' = respSTD;']);
        eval(['A.respAmpE' accelString{i} ' = respAmpE;']);
        eval(['A.respAmpS' accelString{i} ' = respAmpS;']);
    end
end


