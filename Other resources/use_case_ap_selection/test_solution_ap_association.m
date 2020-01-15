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
rng(1899)

load('nn_function.mat');

%Main options
num_it = 50;                         %Number of simulations (random STA positions)
random_pos = 1;                     %Random position of STAs (0:OFF/1:ON)
deploy = 0;                         %Type of deployment (0:circle/1:rectangle)
random_ch = 0;                      %Random selection of the channel
random_traffic = 1;                 %Random traffic selector

%Traffic
L = 12000;                          %Packet length [bits]
lambda_gen = 500;                   %Generation rate per STA [packets/s]
eps_gen = 0;                        %Standard deviation [packets/s]

%Radio module
Pt = 20;                            %Transmission power [dBm]
Sens = -90;                         %Receiver's sensitivity [dBm]
f_backbone = 5E9;                   %Frequency band of backbone links [Hz]
f_access = 2.4E9;                   %Frequency band of access links [Hz]
PL_backbone = 3;                    %Path loss model selector of backbone links
PL_access = 3;                      %Path loss model selector of access links
%PL = 0 (default): Boris model
%PL = 1: ITU-R model
%PL = 2: TMB model (Only for 5 GHz)
%PL = 3: IEEE 802.11ax Residential
%PL = 4: IEEE 802.11ax Enterprise

%IEEE 802.11k/v capabilities (%)
kv_share = [100 0 100 100];
%AP/Extender selection mechanism
score_mode = [99 0 100 101];
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

channel_load_ext = [0 0 0 0 0 0 0 0 0 0 0];

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

if  (random_ch == 0)
    channel_R_a = [1 1 6 11];
end
%AP and extenders parameters and position
d_max_itu = 186.5656;
d_max_11ax_res = 39.7338;
pos_AP = [0 0];
pos_R_a = [0 d_max_11ax_res/4; cos(d_max_11ax_res/4)*d_max_11ax_res/4 ...
    sin(d_max_11ax_res/4)*d_max_11ax_res/4; -cos(d_max_11ax_res/4)*d_max_11ax_res/4 ...
    sin(d_max_11ax_res/4)*d_max_11ax_res/4];

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
M = zeros(dev_a,13);              % WITH Extenders (To choose)

M(1,1) = pos_AP(1);
M(1,2) = pos_AP(2);

M(2:end,1) = pos_R_a(:,1);
M(2:end,2) = pos_R_a(:,2);

% Random distribution of channels
if (random_ch == 1)
    channel_R_a = ChannelDist(f_access,dev_a);
end
% Allocation of channels into M matrix
M(:,8) = transpose(channel_R_a);
    
%-------------------------------------------------------------------------
% PERFORMANCE METRICS VARIABLES
%-------------------------------------------------------------------------
num_approaches = size(score_mode,2);

%OPTION A: RSSI-based
S_T = cell(1,num_approaches);
E_T = cell(1,num_approaches);
E_T_ok = cell(1,num_approaches);
D_avg = cell(1,num_approaches);
D_max = cell(1,num_approaches);
SS_avg = cell(1,num_approaches);
SS_min = cell(1,num_approaches);
assoc_STA = cell(1,num_approaches);
assoc_STA_AP = cell(1,num_approaches);
assoc_STA_E = cell(1,num_approaches);
for i = 1 : num_approaches
    S_T{i} = zeros(num_it,length(lambda_gen));
    E_T{i} = zeros(num_it,length(lambda_gen));
    E_T_ok{i} = zeros(num_it,length(lambda_gen));
    D_avg{i} = zeros(num_it,length(lambda_gen));
    D_max{i} = zeros(num_it,length(lambda_gen));
    SS_avg{i} = zeros(num_it,length(lambda_gen));
    SS_min{i} = zeros(num_it,length(lambda_gen));
    assoc_STA{i} = zeros(num_it,length(lambda_gen));
    assoc_STA_AP{i} = zeros(num_it,length(lambda_gen));
    assoc_STA_E{i} = zeros(num_it,length(lambda_gen));
end

%-------------------------------------------------------------------------
% ALGORITHM EXECUTION
%-------------------------------------------------------------------------
% f = waitbar(0,'','Name','Test AP Association...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

mae = [];
for k = 1 : num_it
    %waitbar(k/num_it,f)
    fprintf('It %d/%d\n', k, num_it);
    % Compute STAs' position
    if (random_pos == 1)
        if (deploy == 0)                %Deployment: circle
            pos_STA = PosGeneratorCircle(pos_AP,sta,R,margin_over,margin_under);
        else                            %Deployment: rectangle
            pos_STA = PosGeneratorRectangle([25 0],sta,side_h,side_v);
        end
    end             
    % Traffic generation ratio (HARDCODED)
    traffic = randi([416,833],sta,1);   
    % N = [posX    posY    Parent_Index    Channel     Rate    Lambda  DBPS     Type    Score_mode]
    %      1       2       3               4           5       6       7        8       9
    N = cell(1, num_approaches);
    for a = 1 : num_approaches    
        %fprintf('- ap %d/%d\n', a, num_approaches);
        N{a} = zeros(sta,9);                 % WITH Extenders (To choose)
        % Allocation of positions into N matrix
        N{a}(:,1) = pos_STA(:,1);
        N{a}(:,2) = pos_STA(:,2);
        % Traffic generation
        N{a}(:,6) = traffic;            
        %STA score mode
        [N{a}(:,8),N{a}(:,9)] = kvGenerator(sta,kv_share(a),score_mode(a));
        % WIFIX Computing framework for AP selection with each of the approaches
        [S_T{a},E_T{a}(k),E_T_ok{a}(k),D_avg{a}(k),D_max{a}(k),SS_avg{a}(k),SS_min{a}(k),assoc_STA{a}(k),assoc_STA_AP{a}(k),assoc_STA_E{a}(k),NEW_MATRIX] = ...
            WIFIXComputing(net, map_R,map_STA,M,N{a},L,Pt,Sens,f_backbone,f_access,PL_backbone,PL_access,ext_conn_alg,score_mode(a),w_a_a,w_b_a,w_c_a,channel_load_ext);
        
        satisfaction_ratio{a}(k,:) = NEW_MATRIX(:,8) ./ N{a}(:,6);
        
        mean_throughput_stas(a, k) = mean(NEW_MATRIX(:,8).*L);
        std_throughput_stas(a, k) = std(NEW_MATRIX(:,8).*L);
        
        if a == 4
        
            mae = [mae; NEW_MATRIX(:,9) * L * 1e-6];
            
        end
        
    end
end

satisf_rand = satisfaction_ratio{1};
satisf_rssi = satisfaction_ratio{2};
satisf_lr = satisfaction_ratio{3};
satisf_nn = satisfaction_ratio{4};


%%
figure
title('Sopoto')
x = [satisf_rand(:) satisf_rssi(:) satisf_lr(:) satisf_nn(:)];
boxplot(x)
% delete(f)

%% PLOTS
figure
mean_throughput = mean(mean_throughput_stas,2)';
std_throughput = mean(std_throughput_stas,2)';
bar(1:size(mean_throughput,2),mean_throughput)
hold on
er = errorbar(1:size(mean_throughput,2),mean_throughput, std_throughput, 'x');
xlabel('Approach','interpreter','latex')
ylabel('$\overline{\Gamma}$ [bps]','interpreter','latex')
xticks(1:size(mean_throughput,2)) 
xticklabels({'Rand', 'RSSI', 'LR', 'NN'})
set(gca,'FontSize',18,'TickLabelInterpreter','latex')

% Percentage of times the NN approach improves the others or remains the same
times_NN_improves_random = (sum(mean_throughput_stas(4,:)>=mean_throughput_stas(1,:))/num_it)*100;
times_NN_improves_RSSI = (sum(mean_throughput_stas(4,:)>=mean_throughput_stas(2,:))/num_it)*100;
times_NN_improves_LR = (sum(mean_throughput_stas(4,:)>=mean_throughput_stas(3,:))/num_it)*100;