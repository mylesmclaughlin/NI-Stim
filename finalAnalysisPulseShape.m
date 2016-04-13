function P = finalAnalysisPulseShape

dataPath = {'C:\Users\Public\Data\NIData\Rat-PulseShape-16-12-2015\'...
    'C:\Users\Public\Data\NIData\Rat-PulseShape-17-12-2015\'...
    'C:\Users\Public\Data\NIData\Rat-Pulse-Shape-10-02-2016\'...
    'C:\Users\Public\Data\NIData\Rat-Pulse-Shape-11-02-2016\'...
    'C:\Users\Public\Data\NIData\Rat-Pulse-Shape-12-02-2016\'};

fileListPS1 = {'Rat-PulseShape-16-Dec-2015-12-07-19.NISmacro'...
    'Rat-PulseShape-17-Dec-2015-11-01-03.NISmacro'...
    'Rat-PulseShape-10-Feb-2016-10-08-40.NISmacro'...
    'Rat-PulseShape-11-Feb-2016-10-30-46.NISmacro'...
    'Rat-PulseShape-12-Feb-2016-09-53-00.NISmacro'};

fileListPS2 = {'Rat-PulseShape-16-Dec-2015-12-22-37.NISmacro'...
    'Rat-PulseShape-17-Dec-2015-12-02-51.NISmacro'...
    'Rat-PulseShape-10-Feb-2016-11-00-31.NISmacro'...
    'Rat-PulseShape-11-Feb-2016-11-26-57.NISmacro'...
    'Rat-PulseShape-12-Feb-2016-10-51-48.NISmacro'};


fileListIPG1 = {''...
    'Rat-IPG-17-Dec-2015-11-21-56.NISmacro'...
    'Rat-PMIPG-10-Feb-2016-10-41-35.NISmacro'...
    'Rat-IPG-11-Feb-2016-10-48-24.NISmacro'...
    'Rat-IPG-12-Feb-2016-10-11-38.NISmacro'};

fileListIPG2 = {''...
    'Rat-IPG-17-Dec-2015-12-20-04.NISmacro'...  %Responses seem smaller on last 2 reps
    ''...
    'Rat-IPG-11-Feb-2016-11-47-49.NISmacro'...
    'Rat-IPG-12-Feb-2016-11-14-06.NISmacro'};

fileListPM1 = {''...
    'Rat-PMIPG-17-Dec-2015-12-40-42.NISmacro'...
    'Rat-PMIPG-10-Feb-2016-11-18-43.NISmacro'...
    'Rat-IPG-11-Feb-2016-11-07-37.NISmacro'...
    'Rat-IPG-12-Feb-2016-11-31-45.NISmacro'};

fileListPM2 = {''...
    ''...
    ''...
    'Rat-PMIPG-11-Feb-2016-12-06-05.NISmacro'...
    ''};

for n = 1:5
    [A,D] = NImakeMacroAvg(fileListPS1{n},dataPath{n});
    P(n) = NIplotMacro(A);
end



figure
list =  {'BP','IPG','PM','TRI','GAUS'};
copo = {'b','g','r','m','y'};
sypo = {'.','.','.','.','.'};
for n = 1:5
    bp_data = squeeze(P(n).respAll(1,:,:));
    ipg_data = squeeze(P(n).respAll(2,:,:));
    pm_data = squeeze(P(n).respAll(3,:,:));
    tri_data = squeeze(P(n).respAll(4,:,:));
    gaus_data = squeeze(P(n).respAll(5,:,:));
    
    allData = [bp_data(:); ipg_data(:); pm_data(:); tri_data(:); gaus_data(:)];
    data_min = min(allData);
    data_max = max(allData-data_min)
    
    bp_data = (bp_data-data_min)/data_max;
    ipg_data = (ipg_data-data_min)/data_max;
    pm_data = (pm_data-data_min)/data_max;
    tri_data = (tri_data-data_min)/data_max;
    gaus_data = (gaus_data-data_min)/data_max;

    subplot(2,2,1)
    plot(bp_data,ipg_data,[copo{n} sypo{n}])
    hold on
    
    subplot(2,2,2)
    plot(bp_data,pm_data,[copo{n} sypo{n}])
    hold on
    
    subplot(2,2,3)
    plot(bp_data,tri_data,[copo{n} sypo{n}])
    hold on
    
    subplot(2,2,4)
    plot(bp_data,gaus_data,[copo{n} sypo{n}])
    hold on
    
    
end
data_min = 0;
data_max = 1;

subplot(2,2,1)
plot([data_min data_max],[data_min data_max],'k:')
xlabel('BP')
ylabel('IPG')
xlim([data_min data_max])
ylim([data_min data_max])

subplot(2,2,2)
plot([data_min data_max],[data_min data_max],'k:')
xlabel('BP')
ylabel('PM')
xlim([data_min data_max])
ylim([data_min data_max])

subplot(2,2,3)
plot([data_min data_max],[data_min data_max],'k:')
xlabel('BP')
ylabel('TRI')
xlim([data_min data_max])
ylim([data_min data_max])

subplot(2,2,4)
plot([data_min data_max],[data_min data_max],'k:')
xlabel('BP')
ylabel('GAUS')
xlim([data_min data_max])
ylim([data_min data_max])

