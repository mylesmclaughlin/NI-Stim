% This functions converts the data of a specific datafile from volts to g
% and acceleration

% This function needs the output 'Callibration of the function
% 'callibration.m' and the output structure S of 'ReadData.m' as the
% input. 

% This function gives the raw data in m/s² and the acceleration in function
% of g as output


function [RawData,Accel_g] = Conversion(Callibration, data)
for j = 1:length(data.dataAccelRaw)
       
        Accel_g.data(j,1) = (data.dataAccelRaw(j,1)-Callibration.x_0g)/Callibration.g_x;
        Accel_g.data(j,2) = (data.dataAccelRaw(j,2)-Callibration.y_0g)/Callibration.g_y;
        Accel_g.data(j,3) = (data.dataAccelRaw(j,3)-Callibration.z_0g)/Callibration.g_z;
        RawData.Time = data.time;
        
        RawData.Acceleration(j,3) = Accel_g.data(j,3)*9.81;
        RawData.Acceleration(j,1) = Accel_g.data(j,1)*9.81;
        RawData.Acceleration(j,2) = Accel_g.data(j,2)*9.81;
end


