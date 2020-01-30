%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% mainS.m --> Main program to test a single network configuration
%-------------------------------------------------------------------------
% INPUT PARAMETERS
%-------------------------------------------------------------------------
%Main options
with_extenders = 1;                 %With extenders (0:OFF/1:ON)
without_extenders = 1;              %Without extenders (0:OFF/1:ON)
random_pos = 1;                     %Random position of STAs (0:OFF/1:ON)
deploy = 0;                         %Type of deployment (0:circle/1:rectangle)
random_traffic = 0;                 %Random traffic selector

%Traffic
L = 12000;                          %Packet length [bits]
lambda_gen = 300;                   %Generation rate per STA [packets/s]

%IEEE 802.11k/v capabilities (%)
kv_share = 100;                     %Share of STAs with 11k/v capabilities (%)

%AP/Extender selection mechanism
score_mode = 5;                     %Score mode employed
%0: 802.11 RSSI-based
%1: Hops & children
%2: Optimized (weights)
%3: Threshold-based (thresholds)
%4: Load aware (paper)
%5: End-to-end latency based (paper)

%OPTION A: WITH Extenders (To choose)
if (score_mode == 0)
    w_a = 0;                        %Not used
    w_b = 0;                        %Not used
    w_c = 0;                        %Not used
    kv_share = 0;
elseif (score_mode == 1)
    w_a = 0;                        %Not used
    w_b = 2;                        %Penalization per hop
    w_c = 4;                        %Penalization per child
elseif (score_mode == 2)
    w_a = 3;                        %Penalization per hop
    w_b = 1000;                     %Penalization per load in the access link
    w_c = 250;                      %Penalization per load in the backbone link
elseif (score_mode == 3)
    w_a = 3;                        %Hop penalization (dB)
    w_b = 70;                       %Maximum load in the access link (%)
    w_c = 50;                       %Maximum load in the backbone link (%)
elseif (score_mode == 4)
    w_a = 0;                        %Not used
    w_b = 1;                        %Penalization per load in the access link
    w_c = 1;                        %Penalization per load in the backbone link
elseif (score_mode == 5)
    w_a = -75;                      %Threshold to directly select the AP
    w_b = 0;                        %Not used
    w_c = 0;                        %Not used
else
    disp('Score mode of OPTION A was not properly chosen');
    kv_share = 0;
end

%Extender connection algorithm
ext_conn_alg = 0;                   %Extender connection algorithm
%ext_conn_alg = 0: Based on RSSI (FON method)
%ext_conn_alg = 1: Based on total delay (UPF method)
margin_R = -70;                     %Margin for direct connection to the AP (dBm)

%Operating Frequency % Channel
WIFI_std = 0;                       %IEEE 802.11 standards employed
%WIFI_std = 0: IEEE 802.11n (2.4 GHz) + IEEE 802.11ac (5 GHz)
%WIFI_std = 1: IEEE 802.11ax (2.4 GHz) + IEEE 802.11ax (5 GHz)
f_backbone = 5E9;                   %Frequency band of backbone links [Hz]
f_access = 2.4E9;                   %Frequency band of access links [Hz]
PL_backbone = 1;                    %Path loss model selector of backbone links
PL_access = 1;                      %Path loss model selector of access links
%PL = 0 (default): Boris model (¿?)
%PL = 1: ITU-R model
%PL = 2: TMB model (Only for 5 GHz)
%PL = 3: IEEE 802.11ax Residential
%PL = 4: IEEE 802.11ax Enterprise

channel_R = [1 6 11 6 11];

%Radio module
Pt = 20;                            %Transmission power [dBm]
Sens = -90;                         %Receiver's sensitivity [dBm]

%STA parameters and position
pos_AP = [0 0];
pos_R = [50 0; 0 -50; -50 0; 0 50];
max_R_per_R = 10;                   %Maximum number of Extenders per Extender
max_STA_per_R = 30;                 %Maximum number of STAs per Extender

if (random_pos == 1)
    sta = 10;                       %Number of stations
    if (deploy == 0)                %Deployment: circle
        margin_under = 0;           %Position margin to deploy STAs far from coverage center [%]
        margin_over = 0;            %Position margin to deploy STAs far from coverage edge [%]
        R = MaxDistance(Pt-Sens,f_access,PL_access);
        pos_STA = PosGeneratorCircle(pos_AP,sta,R,margin_over,margin_under);
    else                            %Deployment: rectangle
        side_h = 100;
        side_v = 20;
        pos_STA = PosGeneratorRectangle(pos_AP,sta,side_h,side_v);
    end
else
    pos_STA = [15 40; -60 40; -10 70; 40 -50; -30 -30; 50 45; -60 -10; -20 -60; -40 45; 70 -25];
    sta = size(pos_STA,1);          %Number of stations
end

%Sizes
exts = size(pos_R,1);               %Number of extenders
dev = exts + 1;                     %Total number of APs + extenders

%Allocating traffic
traffic = zeros(1,sta);             %Traffic pattern [packets/s]
if (random_traffic == 0)
    for i=1:sta
        traffic(i) = lambda_gen;
    end
else
    traffic = [300; 200; 300; 400; 200; 200; 300; 400; 400; 300];
end

% Time Parameters
TPHY = 40E-6;                       %PHY time
SIFS = 16E-6;                       %SIFS time
DIFS = 34E-6;                       %DIFS time
Tslot = 9E-6;                       %Slot time

%Visualization options
map_R = 1;                          %Visualization of Extender map
%0: Do not show Extender map
%1: Show Extender map
map_STA = 0;                        %Visualization of STA map
%0: Do not show STA map
%1: Show final STA map
%2: Show intermediate and final STA maps

%-------------------------------------------------------------------------
% AUXILIAR VARIABLES
%-------------------------------------------------------------------------
% (Note that the first row of M corresponds to the AP, filled with 0s)
% M = [posX    posY    #hops   #Children_R  Parent_Index    Delay   Rate   Channel     #Children_S   lambda_R   DBPS    Airtime access  Airtime backbone]
%      1       2       3       4            5               6       7      8           9             10         11      12              13
M = zeros(dev,13);                  % With Extenders
M_b = zeros(1,13);                  % Without Extenders

% N = [posX    posY    Parent_Index    Channel     Rate    Lambda  DBPS     Type    Score_mode]
%      1       2       3               4           5       6       7        8       9
N = zeros(sta,9);                   % With Extenders
N_b = zeros(sta,9);                 % Without Extenders

routing_table = zeros(dev,1 + max_STA_per_R);

S_STA = zeros(1,sta);
D_STA = zeros(1,sta);
U_STA = zeros(1,sta);
A_STA = zeros(1,sta);

S_R = zeros(1,dev);
D_R = zeros(1,dev);
U_R = zeros(1,dev);
A_R = zeros(1,dev);
SS_R = ones(1,dev);

S_STA_b = zeros(1,sta);
D_STA_b = zeros(1,sta);
U_STA_b = zeros(1,sta);
A_STA_b = zeros(1,sta);

S_T = 0;
S_T_b = 0;

%-------------------------------------------------------------------------
% FILLING MAIN MATRIX
%-------------------------------------------------------------------------
% M :       AP/Extender Matrix (with Extenders)
% M_b:      AP Matrix (without Extenders)
M(1,1) = pos_AP(1);
M(1,2) = pos_AP(2);
M(:,8) = transpose(channel_R);
M_b(1,:) = M(1,:);

M(2:end,1) = pos_R(:,1);
M(2:end,2) = pos_R(:,2);

% N:    STA Matrix
N(:,1) = pos_STA(:,1);
N(:,2) = pos_STA(:,2);
N(:,6) = traffic;

for i=1:sta
    if (rand()*100 < kv_share)
        N(i,8) = 1;
        N(i,9) = score_mode;
    end
end

N_b = N;
M

%% OPTION A: WITH Extenders

close all

if (with_extenders == 1)
    
    disp('*********************************');
    disp('[WITH EXTENDERS] FIRST STAGE: Extenders');
    
    [M] = TopologyExtenders(map_R,M,f_backbone,PL_backbone,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,ext_conn_alg,margin_R,max_R_per_R);
    
    disp('*********************************');
    disp('[WITH EXTENDERS] SECOND STAGE: Stations');
    
    [M,N,routing_table,S_STA,D_STA,U_STA,A_STA,S_R,D_R,U_R,A_R] = TopologySTAs(map_STA,M,N,WIFI_std,f_backbone,f_access,PL_access,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,score_mode,max_STA_per_R,w_a,w_b,w_c);
    
    M
    N
    routing_table
    
    %Auxiliar variables
    active_STA = find(N(:,5) ~= 0);
    direct_STA = transpose(find((N(:,3) == 1)));
    R_hop1 = find((M(:,3) == 1));
    
    %Performance metrics
    S_T = sum(S_STA(direct_STA)) + sum(S_R(R_hop1));        %Total throughput
    E_T = S_T * 100/(sum(traffic) * L);                     %Throughput efficiency
    D_avg = mean(D_STA(active_STA));                        %Average Delay
    SS_STA = S_STA./(transpose(N(:,6))*L)*100;              %Satisfaction [STA]
    SS_R(2:dev) = S_R(2:dev)./(transpose(M(2:dev,10))*L);   %Satisfaction [Extender]
    SS_R = SS_R.*100;
    assoc_STA = sum(N(:,3) ~= 0);                           %Associated STAs
end

%% OPTION B: WITHOUT Extenders
if (without_extenders == 1)
    
    disp('*********************************');
    disp('[WITHOUT EXTENDERS]');
    
    [M_b,N_b,routing_table_b,S_STA_b,D_STA_b,U_STA_b,A_STA_b] = TopologySTAs(map_STA,M_b,N_b,WIFI_std,f_backbone,f_access,PL_access,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,score_mode,max_STA_per_R,0,0,0);
    
    M_b
    N_b
    
    %Auxiliar variables
    active_STA_b = find(N_b(:,5) ~= 0);
    
    %Performance metrics
    S_T_b = sum(S_STA_b(active_STA_b));                 %Total throughput
    E_T_b = S_T_b * 100/(sum(traffic) * L);             %Throughput efficiency
    D_avg_b = mean(D_STA_b(active_STA_b));              %Average Delay
    SS_STA_b = S_STA_b./(transpose(N_b(:,6))*L)*100;    %Satisfaction [STA]
    assoc_STA_b = sum(N_b(:,3) ~= 0);                   %Associated STAs
end

%% RESULTS SUMMARY
disp('*********************************');
disp('[TOTAL THROUGHPUT AT THE AP]');
S_T
S_T_b
disp('THROUGHPUT EFFICIENCY');
E_T
E_T_b
% disp('[STA AVERAGE SATISFACTION]');
% SS_mean = mean(SS_STA)
% SS_mean_b = mean(SS_STA_b)
% disp('[STA MINIMUM SATISFACTION]');
% SS_min = min(SS_STA)
% SS_min_b = min(SS_STA_b)
disp('[STA AVERAGE DELAY]');
D_avg
D_avg_b
disp('[STA MAXIMUM DELAY]');
D_max = max(D_STA)
D_max_b = max(D_STA_b)
disp('[ASSOCIATED STAS]');
assoc_STA
assoc_STA_b
disp('[ASSOCIATED STAS (%)]');
(assoc_STA/sta)*100
(assoc_STA_b/sta)*100

%% COMBINED PLOTS
if ((with_extenders == 1) && (without_extenders == 1))
    
    for i=1:sta
        label_test(i)={[num2str(i),' (CH#',num2str(N(i,4)),')']};
    end
    
    %Maximum rate available [STA]
    Rate_fig = zeros(sta,2);
    Rate_fig(:,1) = N(:,5);
    Rate_fig(:,2) = N_b(:,5);
    
    figure
    bar(Rate_fig/1e6);
    title('Maximum rate available [Mbps]');
    xlabel('STAs')
    xticklabels(label_test)
    legend('WITH Extenders','WITHOUT Extenders');
    
    %Throughput [STA]
    S_fig = zeros(sta,2);
    S_fig(:,1) = S_STA;
    S_fig(:,2) = S_STA_b;
    
    figure
    bar(S_fig/1e6);
    title('Throughput [Mbps]');
    xlabel('STAs')
    xticklabels(label_test)
    legend('WITH Extenders','WITHOUT Extenders');
    
    %Total Throughput [AP]
    S_T_fig = zeros(1,2);
    S_T_fig(1,1) = S_T;
    S_T_fig(1,2) = S_T_b;
    
    figure
    b = bar(S_T_fig/1e6);
    title('Total Throughput in the AP [Mbps]');
    xticklabels({'WITH Extenders','WITHOUT Extenders'})
    
    b.FaceColor = 'flat';
    b.CData(2,:) = [217/255 83/255 25/255];
    
    %Throughput Efficiency [AP]
    E_T_fig = zeros(1,2);
    E_T_fig(1,1) = E_T;
    E_T_fig(1,2) = E_T_b;
    
    figure
    b = bar(E_T_fig);
    title('Throughput efficiency [%]');
    xticklabels({'WITH Extenders','WITHOUT Extenders'})
    
    b.FaceColor = 'flat';
    b.CData(2,:) = [217/255 83/255 25/255];
    
    %     %Satisfaction [STA]
    %     SS_fig = zeros(sta,2);
    %     SS_fig(:,1) = SS_STA;
    %     SS_fig(:,2) = SS_STA_b;
    %
    %     figure
    %     bar(SS_fig);
    %     title('Satisfaction [%]');
    %     xlabel('STAs')
    %     xticklabels(label_test)
    %     legend('WITH Extenders','WITHOUT Extenders');
    %
    %     axis = gca;
    %     axis.Position(3) = 0.6;
    %     annotation('textbox', [0.75, 0.15, 0.1, 0.1], 'String', "Minimum Satisfaction WITH Extenders = " + min(SS_STA) + "%");
    %     annotation('textbox', [0.75, 0.1, 0.1, 0.1], 'String', "Minimum Satisfaction WITHOUT Extenders = " + min(SS_STA_b) + "%");
    %     annotation('textbox', [0.75, 0.05, 0.1, 0.1], 'String', "Average Satisfaction WITH Extenders = " + mean(SS_STA) + "%");
    %     annotation('textbox', [0.75, 0, 0.1, 0.1], 'String', "Average Satisfaction WITHOUT Extenders = " + mean(SS_STA_b) + "%");
    %
    
    %Delay [STA]
    D_fig = zeros(sta,2);
    D_fig(:,1) = D_STA;
    D_fig(:,2) = D_STA_b;
    
    figure
    bar(D_fig);
    title('Total Delay [s]');
    xlabel('STAs')
    xticklabels(label_test)
    legend('WITH Extenders','WITHOUT Extenders');
    
    axis = gca;
    axis.Position(3) = 0.6;
    annotation('textbox', [0.75, 0.15, 0.1, 0.1], 'String', "Maximum Delay WITH Extenders = " + max(D_STA)*1000 + "ms");
    annotation('textbox', [0.75, 0.1, 0.1, 0.1], 'String', "Maximum Delay WITHOUT Extenders = " + max(D_STA_b)*1000 + "ms");
    annotation('textbox', [0.75, 0.05, 0.1, 0.1], 'String', "Average Delay WITH Extenders = " + D_avg*1000 + "ms");
    annotation('textbox', [0.75, 0, 0.1, 0.1], 'String', "Average Delay WITHOUT Extenders = " + D_avg_b*1000 + "ms");
    
    %Utilization [STA]
    U_fig = zeros(sta,2);
    U_fig(:,1) = U_STA;
    U_fig(:,2) = U_STA_b;
    
    figure
    bar(U_fig);
    title('Utilization [%]');
    xlabel('STAs')
    xticklabels(label_test)
    legend('WITH Extenders','WITHOUT Extenders');
    
    %Airtime [STA]
    A_fig = zeros(sta,2);
    A_fig(:,1) = A_STA;
    A_fig(:,2) = A_STA_b;
    
    figure
    bar(A_fig);
    title('Airtime [%]');
    xlabel('STAs')
    xticklabels(label_test)
    legend('WITH Extenders','WITHOUT Extenders');
    
    if (with_extenders == 1)
        
        %Maximum rate available [Extender]
        figure
        bar(M(:,7)/1e6);
        title('Maximum rate available [Mbps]');
        xlabel('AP/Extenders')
        
        %Cumulated traffic [Extender]
        %%%From STAs and possible 'children' extenders
        figure
        bar(M(:,10)*L/1e6);
        title('Cumulated traffic [Mbps]');
        xlabel('AP/Extenders')
        
        %Throughput [Extender]
        S_fig2 = transpose(S_R);
        
        figure
        bar(S_fig2/1e6);
        title('Throughput [Mbps]');
        xlabel('AP/Extenders')
        
        %Satisfaction [Extender]
        SS_fig2 = transpose(SS_R);
        
        figure
        bar(SS_fig2);
        title('Satisfaction [%]');
        xlabel('AP/Extenders')
        
        %Delay [Extender]
        D_fig2 = transpose(D_R);
        
        figure
        bar(D_fig2);
        title('Delay [s]');
        xlabel('AP/Extenders')
        
        %Utilization [Extender]
        U_fig2 = transpose(U_R);
        
        figure
        bar(U_fig2);
        title('Utilization [%]');
        xlabel('AP/Extenders')
        
        %Airtime [Extender]
        A_fig2 = transpose(A_R);
        
        figure
        bar(A_fig2);
        title('Airtime [%]');
        xlabel('AP/Extenders')
        
    end
end

