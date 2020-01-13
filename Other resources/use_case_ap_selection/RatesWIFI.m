%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% RatesWIFI.m --> Computation of data bits per symbol (DBPS) and rate of a
%                 WiFi communication based on received power Pr,
%                 sensitivity S and frequency f
% 2.4 GHz   --> IEEE 802.11n
% 5 GHz     --> IEEE 802.11ac
%-------------------------------------------------------------------------

function [DBPS,rate] = RatesWIFI(Pr,Sens,f)

SD = 52;            %Data Subcarriers
SS = 2;             %Assume 2 Spatial Streams
BPS = 0;            %Number of coded bits per signal carrier for each spatial stream
CR = 0;             %coding rate
rate = 0;

%MCS 0
if(Pr >= Sens && Pr < -79)
    BPS = 1;
    CR = 1/2;
    rate = 7.2E6;
end

%MCS 1
if(Pr >= -79 && Pr < -77)
    BPS = 2;
    CR = 1/2;
    rate = 14.4E6;
end

%MCS 2
if(Pr >= -77 && Pr < -74)
    BPS = 2;
    CR = 3/4;
    rate = 21.7E6;
end

%MCS 3
if(Pr >= -74 && Pr < -70)
    BPS = 4;
    CR = 1/2;
    rate = 28.9E6;
end

%MCS 4
if(Pr >= -70 && Pr < -66)
    BPS = 4;
    CR = 3/4;
    rate = 43.3E6;
end

%MCS 5
if(Pr >= -66 && Pr < -65)
    BPS = 6;
    CR = 2/3;
    rate = 57.8E6;
end

%MCS 6
if(Pr >= -65 && Pr < -64)
    BPS = 6;
    CR = 3/4;
    rate = 65E6;
end

%MCS 7
%Només 2.4GHz
if((Pr >= -64) && (f == 2.4E9))
    BPS = 6;
    CR = 5/6;
    rate = 72.2E6;
end

%MCS 7
%Només 5GHz
if((Pr >= -64 && Pr < -59) && (f == 5E9))
    BPS = 6;
    CR = 5/6;
    rate = 72.2E6;
end

%MCS 8
%Només 5GHz
if((Pr >= -59) && (f == 5E9))
    BPS = 8;
    CR = 3/4;
    rate = 86.7E6;
end

DBPS = SD * SS * BPS * CR;
rate = rate * SS;

end
