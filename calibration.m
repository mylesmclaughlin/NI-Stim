% This function callibrates the accelerometer. The accelerometer needs to
% be hold still for 5 seconds in each direction. Z_still, X_still and
% Y_still are the datasets of these 5 seconds calibration. If you don't
% want to use the 2nd to 7nd second of data measured, the 5 seconds frame can
% be further specified via callibration_z, callibration_x, callibration_y.
% In the format [first datapoint, firstdatapoint-1+5*fs]
% Either all or none of the frames must be further specified.

% Calibration data files should be stored in the folder NIData located in
% the MATLAB path. The MATLAB folder should also contain the functions
% 'ReadData.m', 'binappend.m', 'binfileinfo.m', 'binread.m',
% 'binwrite.m'

% As input this function needs the filename without the full path name (eg.
% 'HumanAccel_Zaxis_still.bin') for all three axes in the order: X_still,
% Y_still, Z_still. Optionally the calibration frame can be further
% specified in the same order.

% As output, two structures are given. The first structure 'Calibration'
% contains the values in volts for 1 g and 0 g for the three axes. The
% second structure, 'O1' contains the calibrated calibration data. This can
% be used as a validation.

function [Calibration] = calibration(X_still, Y_still, Z_still, varargin)

% Read datafiles
fileName1 = [cd '\NIData\' Z_still];
S1 = ReadData(fileName1);
fileName2 = [cd '\NIData\' X_still];
S2 = ReadData(fileName2);
fileName3 = [cd '\NIData\' Y_still];
S3 = ReadData(fileName3);

% Sampling frequency
fs_z = S1.fs;
fs_x = S2.fs;
fs_y = S3.fs;

% Extract calibration data
if nargin <4
    callibr_z = 2*fs_z:7*fs_z;
    callibr_x = 2*fs_x:7*fs_x;
    callibr_y = 2*fs_y:7*fs_y;
else
    callibr_z = varargin{1}(1)*fs_z+1:varargin{1}(2)*fs_z;
    callibr_x = varargin{2}(1)*fs_x+1:varargin{2}(2)*fs_x;
    callibr_y = varargin{3}(1)*fs_y+1:varargin{3}(2)*fs_y;
end

Sz_accelRaw = S1.dataAccelRaw(callibr_z,:);
Sx_accelRaw = S2.dataAccelRaw(callibr_x,:);
Sy_accelRaw = S3.dataAccelRaw(callibr_y,:);

% Create values for 0g and 1g 
x_0g_all = vertcat(Sz_accelRaw(:,1),Sy_accelRaw(:,1));
y_0g_all = vertcat(Sz_accelRaw(:,2), Sx_accelRaw(:,2));
z_0g_all = vertcat(Sx_accelRaw(:,3), Sy_accelRaw(:,3));

x_1g = mean(Sx_accelRaw(:,1));
y_1g = mean(Sy_accelRaw(:,2));
z_1g = mean(Sz_accelRaw(:,3));
Calibration.x_0g = mean(x_0g_all);
Calibration.y_0g = mean(y_0g_all);
Calibration.z_0g = mean(z_0g_all);

Calibration.g_x = abs(x_1g-Calibration.x_0g);
Calibration.g_y = abs(y_1g-Calibration.y_0g);
Calibration.g_z = abs(z_1g-Calibration.z_0g);
