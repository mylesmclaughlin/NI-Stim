% This function callibrates the accelerometer. The accelerometer needs to
% be hold still for 5 seconds in each direction. S1, S2 and S3 are the
% datasets for the calibration. S1 equals the data in the X direction, S2
% in the Y direction and S3 in the Z direction. If you don't want to use
% the 2nd to 7nd second of data measured, the 5 seconds frame can be
% further specified via callibration_z, callibration_x, callibration_y. In
% the format [first datapoint, firstdatapoint-1+5*fs] Either all or none of
% the frames must be further specified.

% As input this function needs the filename without the full path name (eg.
% 'HumanAccel_Zaxis_still.bin') for all three axes in the order: x
% direction, y direction, z direction. Also the sampling frequency needs to
% be further specified. Optionally the calibration frame can be further
% specified in the same order.

% As output, a structure is given. The  structure 'Calibration'
% contains the values in volts for 1 g and 0 g for the three axes. 

function [Calibration] = calibration(S1, S2, S3, fs, varargin)

% Extract calibration data
if nargin <5
    callibr_x = 2*fs:7*fs;
    callibr_y = 2*fs:7*fs;
    callibr_z = 2*fs:7*fs;
else
    callibr_x = varargin{1}(1)*fs+1:varargin{1}(2)*fs;
    callibr_y = varargin{2}(1)*fs+1:varargin{2}(2)*fs;
    callibr_z = varargin{3}(1)*fs+1:varargin{3}(2)*fs;
end

Sx_accelRaw = S1(callibr_x,:);
Sy_accelRaw = S2(callibr_y,:);
Sz_accelRaw = S3(callibr_z,:);

% Create values for 0g and 1g 
x_0g_all = vertcat(Sy_accelRaw(:,1), Sz_accelRaw(:,1));
y_0g_all = vertcat(Sx_accelRaw(:,2), Sz_accelRaw(:,2));
z_0g_all = vertcat(Sx_accelRaw(:,3), Sy_accelRaw(:,3));

x_1g = mean(Sx_accelRaw(:,1));
y_1g = mean(Sy_accelRaw(:,2));
z_1g = mean(Sz_accelRaw(:,3));

% Make Calibration structure
Calibration.x_0g = mean(x_0g_all);
Calibration.y_0g = mean(y_0g_all);
Calibration.z_0g = mean(z_0g_all);

Calibration.g_x = abs(x_1g-Calibration.x_0g);
Calibration.g_y = abs(y_1g-Calibration.y_0g);
Calibration.g_z = abs(z_1g-Calibration.z_0g);
