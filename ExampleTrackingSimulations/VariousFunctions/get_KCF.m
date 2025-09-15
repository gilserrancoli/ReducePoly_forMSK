% This function returns the KCF given as inputs the path to the
% file containing experimental values of KCF
%
% Author: Gil Serrancolí
% Date: 09/07/2021
% 
function KCF=get_KCF(pathKCF)

load(pathKCF,'dataKCF');
KCF_init.time=dataKCF.data(:,contains(dataKCF.colheaders,{'Time'}));
KCF_init.data(:,1)=KCF_init.time;

KCF_init.data(:,2)=dataKCF.data(:,contains(dataKCF.colheaders,{'PM'}))+...
    dataKCF.data(:,contains(dataKCF.colheaders,{'AM'}));
KCF_init.data(:,3)=dataKCF.data(:,contains(dataKCF.colheaders,{'AL'}))+...
    dataKCF.data(:,contains(dataKCF.colheaders,{'PL'}));

% Low-pass filter
order = 4;
cutoff_low = 6;
fs=1/mean(diff(KCF_init.time));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
KCF.time=KCF_init.time;
KCF.data(:,1)=KCF_init.time;
KCF.data(:,2:3) = filtfilt(af,bf,KCF_init.data(:,2:3))*4.4482216153; %data in csv was in lbs, we convert them to N
end