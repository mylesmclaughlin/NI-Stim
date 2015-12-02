%% Random data
close all
tic
addpath('NIData');
addpath('Functions');
fileName1 = 'C:\Users\Public\Data\NIData\Myles-inv-trem_1.bin';
S = ReadData(fileName1);
toc
% Calibration based on still data of the 3 main axes
tic
[Calibration] = calibration('x_calib01-Dec-2015-11-12-06.bin', 'y_calib01-Dec-2015-11-12-23.bin', 'z_calib01-Dec-2015-11-12-40.bin');
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

%% Frequency and phase
% figure
% fftDisp = mfft(Disptot,S.fs);
% xlim([0 10])
% 
% 
% freqz(Disptot)
% figure
% y = hilbert(Disptot);
% plot(Kalman.Time,angle(y))
% 
% figure
% y = hilbert(S.current);
% plot(Kalman.Time,angle(y))
% 
% figure
% fftCur = mfft(S.current,S.fs);
% xlim([0 10]);

stimperiod = 1/9; %stimfreq = 9Hz
stimperiodsamp = round(stimperiod*S.fs);

DisptotNoMean = Disptot-mean(Disptot);
[corr, lags] = xcorr(S.current,DisptotNoMean);

figure;
plot(corr(49400:end,:))
hold on
plot(corr)

figure
plot(lags)

%%
figure
spectrogram(DisptotNoMean,1024,500, 1024, 1000,'yaxis')
