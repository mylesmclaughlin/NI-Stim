% This function integrates the raw acceleration to velocity and
% displacement. Input is the RawData which is the output of the Conversion
% function

function [KalmanIterative] = IterativeKalman(RawData, S)

%% Detrend Acceleration

RawAccelDetrend = detrend(RawData.Acceleration,'constant');
KalmanIterative.AccelDetrend = RawAccelDetrend;
KalmanIterative.Time = RawData.Time;

%% Apply KalmanIterative filter for individual axes

% Displacement
A = [1 1/S.fs; 0 1];
B = [0;1/S.fs];
C = [1 0];
D = [0];
ts = 1/S.fs;
FilterSys = ss(A,B,C,D,ts,'InputName','Acceleration', 'OutputName', 'position');
sampVec = [1:S.fs];
KalmanIterative.DisplacementZ = [];
KalmanIterative.DisplacementX = [];
KalmanIterative.DisplacementY = [];
KalmanIterative.SpeedZ = [];
KalmanIterative.SpeedX = [];
KalmanIterative.SpeedY = [];
Q = [1]; R = [1];
x0x = [0; 0];
x0y = [0; 0];
x0z = [0; 0];
 for i = 1:floor(length(RawAccelDetrend)/S.fs)
    [kalmfdisp, L,P,M] = kalman(FilterSys,Q,R,0);

    [outz, t, z] = lsim(kalmfdisp, [RawAccelDetrend(sampVec,3)], 0:1/S.fs:1-1/S.fs, x0z);
    KalmanIterative.DisplacementZ = [KalmanIterative.DisplacementZ; z(:,1)];
    KalmanIterative.SpeedZ = [KalmanIterative.SpeedZ; z(:,2)]; % Z gebruiken want ik heb puur de states nodig als output
    [outx, t, x] = lsim(kalmfdisp, [RawAccelDetrend(sampVec,1)], 0:1/S.fs:1-1/S.fs, x0x);
    KalmanIterative.DisplacementX = [KalmanIterative.DisplacementX ; x(:,1)];
    KalmanIterative.SpeedX = [KalmanIterative.SpeedX; x(:,2)];
    [outy, t, y] = lsim(kalmfdisp, [RawAccelDetrend(sampVec,2)], 0:1/S.fs:1-1/S.fs, x0y);
    KalmanIterative.DisplacementY = [KalmanIterative.DisplacementY; y(:,1)];
    KalmanIterative.SpeedY = [KalmanIterative.SpeedY; y(:,2)];
    sampVec = sampVec + S.fs;
    x0z = z(end,:);
    x0y = y(end,:);
    x0x = x(end,:);
 end
 
 KalmanIterative.Time = RawData.Time(1,1:(floor(length(RawAccelDetrend)/S.fs))*S.fs);
