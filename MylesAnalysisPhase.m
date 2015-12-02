%% Add data
close all
tic
addpath('NIData');
addpath('Functions');
fileName1 = [cd '\NIData\Myles.bin'];
S = ReadData(fileName1);
toc
% Calibration based on still data of the 3 main axes
tic
[Calibration] = calibration('x_calib24-Nov-2015-16-05-00.bin', 'y_calib24-Nov-2015-16-05-16.bin', 'z_calib24-Nov-2015-16-05-33.bin');
toc

% Calibration data
tic
[RawData,~] = Conversion(Calibration, S);
toc
% Integration to velocity and displacement
tic
[Kalman] = IterativeKalman(RawData,S);
toc


%% Figures

%% Total displacement
Disptot = sqrt(Kalman.DisplacementX.^2+Kalman.DisplacementY.^2+Kalman.DisplacementZ.^2);

figure
plot(Kalman.Time,Kalman.DisplacementX)
hold on 
plot(Kalman.Time,Kalman.DisplacementY)
hold on
plot(Kalman.Time,Kalman.DisplacementZ)
legend('X','Y','Z')

figure
plot(Kalman.Time,Disptot)

figure
plot(S.time,S.current)

%% General phase
stimperiod = 1/9; %stimfreq = 9Hz
stimperiodsamp = round(stimperiod*S.fs); % aantal samples in 1 periode van stimulatie

DisptotNoMean = Disptot-mean(Disptot);
[corr, lags] = xcorr(S.current,DisptotNoMean,stimperiodsamp);
figure; plot(lags/S.fs,corr); %phase difference biggest around 0.49s after stimulation
xlabel('Time after current stimulation');ylabel('Crosscorrelation');axis tight

%% Phase over time

CorrIterative = [];
lagsIterative = [];
for i = 1:length(Disptot)-1*S.fs
    CorrData = DisptotNoMean(i:i+1*S.fs);
    CorrCurrent = S.current(i:i+1*S.fs);
    [corr,~] = xcorr(CorrCurrent,CorrData,stimperiodsamp,'unbiased');
    [correlation,lags] = max(abs(corr));
    CorrIterative = [CorrIterative; corr(lags)];
    lagsIterative = [lagsIterative;(lags-stimperiodsamp)/S.fs];
end
figure;
subplot(2,1,1)
plot(1/S.fs:1/S.fs:(length(Disptot/S.fs)-1*S.fs)/S.fs,lagsIterative); 
xlabel('Time');ylabel('Phase difference (s)');axis tight
title('Phase difference over time')

subplot(2,1,2);
plot(1/S.fs:1/S.fs:(length(Disptot/S.fs)-1*S.fs)/S.fs,CorrIterative); 
xlabel('Time');ylabel('Crosscorrelation');axis tight
title('Crosscorrelation over time')