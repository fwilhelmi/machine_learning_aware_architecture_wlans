%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% PacketDelay.m --> Computation of the delay of a single WiFi packet
%-------------------------------------------------------------------------

function T = PacketDelay(L,TPHY,DBPS,SIFS,DIFS,Tslot)

T = TPHY + ceil((16+224+L+6) ./ DBPS).*4E-6+SIFS+ceil((16+112+L+6) ./ DBPS).*4E-6 + DIFS + 7.5.*Tslot;

end