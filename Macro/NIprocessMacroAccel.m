function D = NIprocessMacroAccel(D,accelInd)

accelString = {'X','Y','Z'};
bandpass = [1 500];
filterorder = 2;
[filterb filtera] = butter(filterorder,bandpass/(A.fs/2),'bandpass');
zInd = find(D.tvec<0);
m = 1;

for a = 1:length(accelInd)    
    for s = 1:D.nSeq
        stimStartInd = zInd(end)+1;
        eval(['allData = squeeze(D.allData' D.macroList{m} '(:,s,:,accelInd(a)));'])
        
        ymax = [];
        ymin = [];
        for r = 1:D.nReps
            calib = mean(allData(r,zInd));
            allData(r,:) = allData(r,:) - calib;
            allData(r,:) = filter(filterb,filtera,allData(r,:));
            allData(r,:) = cumtrapz(allData(r,:));
            allData(r,:) = filter(filterb,filtera,allData(r,:));
            allData(r,:) = cumtrapz(allData(r,:));
        end
        
        ydata = mean(allData);
        ydata_std = std(allData);
        
        ymax(s) = max(ydata);
        ymin(s) = min(ydata);
        [respAmp(s),ind] = max(abs(ydata(stimStartInd:end)));
        respSTD(s) = ydata_std(ind);
    end
    eval(['A.respAmp' accelString{i} ' = respAmp;']);
    eval(['A.respSTD' accelString{i} ' = respSTD;']);
    eval(['A.respAmpE' accelString{i} ' = respAmpE;']);
    eval(['A.respAmpS' accelString{i} ' = respAmpS;']);
end