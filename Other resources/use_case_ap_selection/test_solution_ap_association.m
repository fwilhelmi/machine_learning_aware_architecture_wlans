%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% main.m --> Main program to test multiple network configurations and
%             plot results
clear
clc
%-------------------------------------------------------------------------
% INPUT PARAMETERS
%-------------------------------------------------------------------------
%Keep the same random generated number
rng(2000)

%Main options
num_it = 1;                      %Number of simulations (random STA positions)
random_pos = 1;                     %Random position of STAs (0:OFF/1:ON)
deploy = 0;                         %Type of deployment (0:circle/1:rectangle)
random_ch = 0;                      %Random selection of the channel
random_traffic = 1;                 %Random traffic selector

%Traffic
L = 12000;                          %Packet length [bits]
%lambda_gen = 1:1:300;              %Generation rate per STA [packets/s]
lambda_gen = 500;                   %Generation rate per STA [packets/s]
%lambda_gen = 1:50:501;             %Generation rate per STA [packets/s]
%lambda_gen = 250;                  %Generation rate per STA [packets/s]
eps_gen = 0;                        %Standard deviation [packets/s]

%IEEE 802.11k/v capabilities (%)
kv_share_a = 0;                   %Share of STAs with 11k/v capabilities (%) in Option A
kv_share_b = 100;                   %Share of STAs with 11k/v capabilities (%) in Option A
kv_share_c = 100;                   %Share of STAs with 11k/v capabilities (%) in Option A

%AP/Extender selection mechanism
score_mode_a = 0;                     %Score mode employed for Option A (WITH Extenders - To choose)
score_mode_b = 99;                    %Score mode employed for Option B (WITH Extenders - To choose)
score_mode_c = 100;                   %Score mode employed for Option C (WITH Extenders - To choose)
%0: 802.11 RSSI-based
%1: Hops & children
%2: Optimized (weights)
%3: Threshold-based (thresholds)
%4: Load aware (paper)
%5: End-to-end latency based (paper)
%99: Random selection (only extenders)
%100: Linear regression approach (only extenders)

w_a_a = 0;                      %Not used
w_b_a = 0;                      %Penalization per hop
w_c_a = 0;                      %Penalization per child

% Time Parameters
TPHY = 40E-6;                       %PHY time
SIFS = 16E-6;                       %SIFS time
DIFS = 34E-6;                       %DIFS time
Tslot = 9E-6;                       %Slot time

channel_load_ext = [0 0 0 0 0 0 0 0 0 0 0];
%Radio module
Pt = 20;                            %Transmission power [dBm]
Sens = -90;                         %Receiver's sensitivity [dBm]

%Computation options
map_R = 0;                          %Visualization of Extender map
%0: Do not show Extender map
%1: Show Extender map
map_STA = 0;                        %Visualization of STA map
%0: Do not show STA map
%1: Show final STA map (it affects computation time!!!)
%2: Show intermediate and final STA maps (it affects computation time!!!)

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
PL_backbone = 3;                    %Path loss model selector of backbone links
PL_access = 3;                      %Path loss model selector of access links
%PL = 0 (default): Boris model (�?)
%PL = 1: ITU-R model
%PL = 2: TMB model (Only for 5 GHz)
%PL = 3: IEEE 802.11ax Residential
%PL = 4: IEEE 802.11ax Enterprise
if  (random_ch == 0)
    channel_R_a = [1 1 6 11];
end
%STA parameters and position
d_max_itu = 186.5656;
d_max_11ax_res = 39.7338;
pos_AP = [0 0];
pos_R_a = [0 d_max_11ax_res/4; cos(d_max_11ax_res/4)*d_max_11ax_res/4 ...
    sin(d_max_11ax_res/4)*d_max_11ax_res/4; -cos(d_max_11ax_res/4)*d_max_11ax_res/4 ...
    sin(d_max_11ax_res/4)*d_max_11ax_res/4];
max_R_per_R = 30;                   %Maximum number of Extenders per Extender
max_STA_per_R = 30;                 %Maximum number of STAs per Extender
if (random_pos == 1)
    sta = 20;                       %Number of stations
    %Deployment: circle
    margin_under = 0;               %Position margin to deploy STAs far from coverage center [%]
    margin_over = 0;                %Position margin to deploy STAs far from coverage edge [%]
    R = MaxDistance(Pt-Sens,f_access,PL_access)/2;
    %Deployment: rectangle
    side_h = 100;
    side_v = 20;
else
    pos_STA = [75 10; -25 -10; 65 5; 30 5; 30 -5; 55 -5; 60 -5; 65 -5; 70 -5; -25 10];
    sta = size(pos_STA,1);          %Number of stations
    num_it = 1;
end
%Sizes
exts_a = size(pos_R_a,1);           %Number of extenders
dev_a = exts_a + 1;                 %Total number of APs + extenders

%-------------------------------------------------------------------------
% AUXILIAR VARIABLES
%-------------------------------------------------------------------------
% (Note that the first row of M corresponds to the AP, filled with 0s)
% M = [posX    posY    #hops   #Children_R  Parent_Index    Delay   Rate   Channel     #Children_S   lambda_R   DBPS    Airtime access  Airtime backbone]
%      1       2       3       4            5               6       7      8           9             10         11      12              13
M_a = zeros(dev_a,13);              % WITH Extenders (To choose)

M_a(1,1) = pos_AP(1);
M_a(1,2) = pos_AP(2);

M_a(2:end,1) = pos_R_a(:,1);
M_a(2:end,2) = pos_R_a(:,2);

f = waitbar(0,'','Name','Test AP Association...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

total_repetitions = 1000;
for r = 1 : total_repetitions
    %-------------------------------------------------------------------------
    % PERFORMANCE METRICS VARIABLES
    %-------------------------------------------------------------------------
    %Definition of variables
    %OPTION A: RSSI-based
    S_T_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    E_T_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    E_T_ok_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    D_avg_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    D_max_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    SS_avg_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    SS_min_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_AP_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_E_a = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));

    S_T_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    E_T_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    share_ok_a = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    D_avg_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    D_max_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    SS_avg_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    SS_min_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_AP_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_E_a_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));

    %OPTION B: Random
    S_T_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    E_T_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    E_T_ok_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    D_avg_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    D_max_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    SS_avg_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    SS_min_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_AP_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_E_b = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));

    S_T_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    E_T_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    share_ok_b = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    D_avg_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    D_max_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    SS_avg_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    SS_min_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_AP_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_E_b_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));

    %OPTION C: Linear regression
    S_T_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    E_T_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    E_T_ok_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    D_avg_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    D_max_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    SS_avg_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    SS_min_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_AP_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));
    assoc_STA_E_c = zeros(length(w_b_a),length(w_c_a),num_it,length(lambda_gen));

    S_T_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    E_T_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    share_ok_c = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    D_avg_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    D_max_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    SS_avg_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    SS_min_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_AP_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));
    assoc_STA_E_c_avg = zeros(length(w_b_a),length(w_c_a),length(lambda_gen));

    % ALGORITHM EXECUTION
    for k = 1:num_it
        % Random distribution of channels
        if (random_ch == 1)
            channel_R_a = ChannelDist(f_access,dev_a);
        end
        % Allocation of channels into M matrix
        M_a(:,8) = transpose(channel_R_a);
        % N = [posX    posY    Parent_Index    Channel     Rate    Lambda  DBPS     Type    Score_mode]
        %      1       2       3               4           5       6       7        8       9
        N_a = zeros(sta,9);                 % WITH Extenders (To choose)
        N_b = zeros(sta,9);                 % WITH Extenders (To choose)
        N_c = zeros(sta,9);                 % WITH Extenders (To choose)
        %STA position
        if (random_pos == 1)
            if (deploy == 0)                %Deployment: circle
                pos_STA = PosGeneratorCircle(pos_AP,sta,R,margin_over,margin_under);
            else                            %Deployment: rectangle
                pos_STA = PosGeneratorRectangle([25 0],sta,side_h,side_v);
            end
        end
        % Allocation of positions into N matrix
        N_a(:,1) = pos_STA(:,1);
        N_a(:,2) = pos_STA(:,2);
        N_b(:,1) = pos_STA(:,1);
        N_b(:,2) = pos_STA(:,2);
        N_c(:,1) = pos_STA(:,1);
        N_c(:,2) = pos_STA(:,2);
        %STA score mode
        [N_a(:,8),N_a(:,9)] = kvGenerator(sta,kv_share_a,score_mode_a);
        [N_b(:,8),N_b(:,9)] = kvGenerator(sta,kv_share_b,score_mode_b);
        [N_c(:,8),N_c(:,9)] = kvGenerator(sta,kv_share_c,score_mode_c);
        %Code loop
        for l = 1:length(lambda_gen)
            m = 0;
            % Traffic model
            if (random_traffic == 0)
                N_a(:,6) = lambda_gen(l);
                N_b(:,6) = lambda_gen(l);
                N_c(:,6) = lambda_gen(l);
            else
                N_a(:,6) = randi([416,833],sta,1);   % HARDCODED
                N_b(:,6) = randi([416,833],sta,1);   % HARDCODED
                N_c(:,6) = randi([416,833],sta,1);   % HARDCODED
    %             for p = 1:sta
    %                 %N_a(p,6) = eps_gen.*randn() + lambda_gen(l);
    %                 N_a(p,6) = randi([1,5],sta,1);
    %             end
            end
            for i = 1:length(w_b_a)
                for j = 1:length(w_c_a)
                    if (score_mode_a == 4)
                        w_c_a(j) = 1 - w_b_a(i);
                    end
                    m = m + 1;
                    waitbar(r/total_repetitions,f)
                    %OPTION A: RSSI-based
                    [S_T_a,E_T_a(i,j,k,l),E_T_ok_a(i,j,k,l),D_avg_a(i,j,k,l),D_max_a(i,j,k,l),SS_avg_a(i,j,k,l),SS_min_a(i,j,k,l),assoc_STA_a(i,j,k,l),assoc_STA_AP_a(i,j,k,l),assoc_STA_E_a(i,j,k,l)] = ...
                        WIFIXComputing(map_R,map_STA,M_a,N_a,WIFI_std,f_backbone,f_access,PL_backbone,PL_access,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,ext_conn_alg,margin_R,max_R_per_R,score_mode_a,max_STA_per_R,w_a_a,w_b_a(i),w_c_a(j),channel_load_ext);
                    %OPTION B: Random
                    [S_T_b,E_T_b(i,j,k,l),E_T_ok_b(i,j,k,l),D_avg_b(i,j,k,l),D_max_b(i,j,k,l),SS_avg_b(i,j,k,l),SS_min_b(i,j,k,l),assoc_STA_b(i,j,k,l),assoc_STA_AP_b(i,j,k,l),assoc_STA_E_b(i,j,k,l)] = ...
                        WIFIXComputing(map_R,map_STA,M_a,N_b,WIFI_std,f_backbone,f_access,PL_backbone,PL_access,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,ext_conn_alg,margin_R,max_R_per_R,score_mode_b,max_STA_per_R,w_a_a,w_b_a(i),w_c_a(j),channel_load_ext);
                    %OPTION C: Linear regression
                    [S_T_c,E_T_c(i,j,k,l),E_T_ok_c(i,j,k,l),D_avg_c(i,j,k,l),D_max_c(i,j,k,l),SS_avg_c(i,j,k,l),SS_min_c(i,j,k,l),assoc_STA_c(i,j,k,l),assoc_STA_AP_c(i,j,k,l),assoc_STA_E_c(i,j,k,l)] = ...
                        WIFIXComputing(map_R,map_STA,M_a,N_c,WIFI_std,f_backbone,f_access,PL_backbone,PL_access,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,ext_conn_alg,margin_R,max_R_per_R,score_mode_c,max_STA_per_R,w_a_a,w_b_a(i),w_c_a(j),channel_load_ext);
                end
            end
        end
    end
%     %Computing averages for the num_rep performed simulations
%     for i = 1:length(w_b_a)
%         for j = 1:length(w_c_a)
%             for l = 1:length(lambda_gen)
%                 %OPTION A: WITH Extenders (To choose)
%                 [S_T_a_avg(i,j,l),E_T_a_avg(i,j,l),share_ok_a(i,j,l),D_avg_a_avg(i,j,l),D_max_a_avg(i,j,l),SS_avg_a_avg(i,j,l),SS_min_a_avg(i,j,l),assoc_STA_a_avg(i,j,l),assoc_STA_AP_a_avg(i,j,l),assoc_STA_E_a_avg(i,j,l)] = ...
%                     MeanGenerator(S_T_a(i,j,:,l),E_T_a(i,j,:,l),E_T_ok_a(i,j,:,l),D_avg_a(i,j,:,l),D_max_a(i,j,:,l),SS_avg_a(i,j,:,l),SS_min_a(i,j,:,l),assoc_STA_a(i,j,:,l),assoc_STA_AP_a(i,j,:,l),assoc_STA_E_a(i,j,:,l),num_it);
%             %OPTION A: WITH Extenders (To choose)
%                 [S_T_b_avg(i,j,l),E_T_b_avg(i,j,l),share_ok_b(i,j,l),D_avg_b_avg(i,j,l),D_max_b_avg(i,j,l),SS_avg_b_avg(i,j,l),SS_min_b_avg(i,j,l),assoc_STA_b_avg(i,j,l),assoc_STA_AP_b_avg(i,j,l),assoc_STA_E_b_avg(i,j,l)] = ...
%                     MeanGenerator(S_T_b(i,j,:,l),E_T_b(i,j,:,l),E_T_ok_b(i,j,:,l),D_avg_b(i,j,:,l),D_max_b(i,j,:,l),SS_avg_b(i,j,:,l),SS_min_b(i,j,:,l),assoc_STA_b(i,j,:,l),assoc_STA_AP_b(i,j,:,l),assoc_STA_E_b(i,j,:,l),num_it);
%             end
%         end
%     end
    solution_a(r) = mean(S_T_a);
    solution_b(r) = mean(S_T_b);
    solution_c(r) = mean(S_T_c);
end

delete(f)

%% PLOTS
figure
mean1 = mean(solution_a);
mean2 = mean(solution_b);
mean3 = mean(solution_c);
std1 = std(solution_a);
std2 = std(solution_b);
std3 = std(solution_c);
bar(1:2,[mean1, mean3])
hold on
er = errorbar(1:2,[mean1, mean3], [std1, std3], 'x');
xlabel('Approach')
ylabel('Mbps')
xticks([1 2]) 
xticklabels({'RSSI-based', 'Linear Regr.'})
set(gca,'FontSize',18)