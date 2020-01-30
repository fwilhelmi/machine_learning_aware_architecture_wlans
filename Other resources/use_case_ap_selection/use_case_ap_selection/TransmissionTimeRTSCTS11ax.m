% Packet Aggregation is not considered
% 11ax

function [Ts,Tc] = TransmissionTimeRTSCTS11ax(Na,L,B,SU,Ym,Yc)

    %Na --> Aggregated packets
    %L -->  Packet length
    %B -->  Bandwidth
    %SU --> Spatial Streams
    %Ym --> Bits per symbol
    %Yc --> Coding rate

    % RTS + SIFS + CTS + SIFS + Data + SIFS + ACK + DIFS + Te;

    EmptySlot = 9E-6;
    SIFS = 16E-6;
    DIFS = 34E-6;

    PHY_SU = 164E-6;
    PHY_Legacy = 20E-6;

    RTS = 160;
    CTS = 112;

    MH = 320;
    SF = 16;
    TB = 18;
    MD = 32;

    BA = 240; % Is the compressed version
    BA = 432; % Per 256;

    % -----------Transmission Rate -------------------
    Ts = 16E-6;
    Ts_legacy = 4E-6;

    Ysb = NumberOfSubcarriers11ax(B);
    BitsOFDM = Ysb*Ym*Yc*SU;
    BitsOFDMLegacyRate = 48*1*(1/2); % 6 Mbps

    % ------------------------------------------------

    % Successful Transmission

    T_RTS = PHY_Legacy + ceil((SF+RTS+TB)/BitsOFDMLegacyRate)*Ts_legacy;
    T_CTS = PHY_Legacy + ceil((SF+CTS+TB)/BitsOFDMLegacyRate)*Ts_legacy;
    T_BA = PHY_Legacy + ceil((SF+BA+TB)/BitsOFDMLegacyRate)*Ts_legacy;
    %T_BA = PHY_SU + ceil((SF+BA+TB)/BitsOFDM)*Ts


    T_DATA = 0;
    NaA=1;
    for na = 1:Na
        if(Na>1) 
            T_DATA_PREV = PHY_SU + ceil((SF+na*(MD+MH+L)+TB)/BitsOFDM)*Ts;
        else
            T_DATA_PREV = PHY_SU + ceil((SF+(MH+L)+TB)/BitsOFDM)*Ts;
        end

        if(T_DATA_PREV > 5.484E-3) 
            %disp(T_DATA_PREV);
            %disp(na); 
            %NaA=na;
            %pause
            if(na==1)
               T_DATA = T_DATA_PREV; 
            end        
            break;
        else
            NaA=na;
            T_DATA = T_DATA_PREV;
        end
    end

    % disp(T_DATA);
    % disp(NaA);
    % pause

    %Ts = T_RTS + SIFS + T_CTS + SIFS + T_DATA + SIFS + T_BAR + SIFS + T_BA + DIFS + Te;

    Ts = T_RTS + SIFS + T_CTS + SIFS + T_DATA + SIFS + T_BA + DIFS + EmptySlot;

    %Ts = T_DATA + SIFS + T_BA + DIFS + Te;

    T_DATAx = (Na*L/BitsOFDM)*Ts; 

    %pause

    % Throughput

    S=Na*L/Ts;

    % Collision

    Tc = T_RTS + SIFS + T_CTS + DIFS + EmptySlot;

    %[Na NaA T]

    %pause
end