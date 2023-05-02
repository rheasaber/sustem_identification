function [y,u,Ts] = GetExperimentData(path_)
if nargin == 0
    path_ = 'logs.bin';
end
FID = fopen(path_, 'r');
DAT = fread(FID,'double');
fclose(FID);


y = DAT(1:6:end);
u = DAT(2:6:end);

flag = logical(DAT(6:6:end)); % flag is '1' when identification signal is applied


y = y(flag);
u = u(flag);

Ts = 10e-3;
end