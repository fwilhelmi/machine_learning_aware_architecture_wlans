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

num_stas = [10 20 30];

%Main options
num_it = 1;                         %Number of simulations (random STA positions)
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
ext_conn_alg = 0;                   % AP connection algorithm
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

sta = 30;                       %Number of stations
%Deployment: circle
margin_under = 0;               %Position margin to deploy STAs far from coverage center [%]
margin_over = 0;                %Position margin to deploy STAs far from coverage edge [%]
R = MaxDistance(Pt-Sens,f_access,PL_access)/2;
%Deployment: rectangle
side_h = 100;
side_v = 20;

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

mae = [];
for k = 1 : num_it
    fprintf('It %d/%d\n', k, num_it);
    % Compute STAs' position
    if (deploy == 0)                %Deployment: circle
        pos_STA = PosGeneratorCircle(pos_R_a(1,:),sta,R,margin_over,margin_under);
    else                            %Deployment: rectangle
        pos_STA = PosGeneratorRectangle([25 0],sta,side_h,side_v);
    end
    DrawStas(pos_STA, pos_R_a);
    % Traffic generation ratio
    traffic = randi([2*83,833],sta,1);   
    % N = [posX    posY    Parent_Index    Channel     Rate    Lambda  DBPS     Type    Score_mode]
    %      1       2       3               4           5       6       7        8       9
    N = cell(1, num_approaches);
    for a = 1 : num_approaches 
        satisfaction_ratio{a} = zeros(num_it, sta);
        N{a} = zeros(sta,9);                 % WITH Extenders (To choose)
        % Allocation of positions into N matrix
        N{a}(:,1) = pos_STA(:,1);
        N{a}(:,2) = pos_STA(:,2);
        % Traffic generation
        N{a}(:,6) = traffic;            
        % STA score mode
        [N{a}(:,8),N{a}(:,9)] = kvGenerator(sta,kv_share(a),score_mode(a));
        % Framework for AP selection with each of the approaches
        [S_T{a},E_T{a}(k),E_T_ok{a}(k),D_avg{a}(k),D_max{a}(k),SS_avg{a}(k),SS_min{a}(k),assoc_STA{a}(k),assoc_STA_AP{a}(k),assoc_STA_E{a}(k),NEW_MATRIX] = ...
            WIFIXComputing(net, map_R,map_STA,M,N{a},L,Pt,Sens,f_backbone,f_access,PL_backbone,PL_access,ext_conn_alg,score_mode(a),w_a_a,w_b_a,w_c_a,channel_load_ext);

        satisfaction_ratio{a}(k,:) = NEW_MATRIX(:,8) ./ N{a}(:,6);

        mean_throughput_stas(a, k) = mean(NEW_MATRIX(:,8).*L);
        std_throughput_stas(a, k) = std(NEW_MATRIX(:,8).*L);

        min_tpt{1}(a,k) = min(S_T{a}(1:8));
        min_tpt{2}(a,k) = min(S_T{a}(9:16));
        min_tpt{3}(a,k) = min(S_T{a}(17:24));
        
        agg_tpt{1}(a,k) = sum(S_T{a}(1:8));
        agg_tpt{2}(a,k) = sum(S_T{a}(9:16));
        agg_tpt{3}(a,k) = sum(S_T{a}(17:24));
        
        tpt{1}(a,k,:) = S_T{a}(1:8);
        tpt{2}(a,k,:) = S_T{a}(9:16);
        tpt{3}(a,k,:) = S_T{a}(17:24);

        if a == 4
            mae = [mae; NEW_MATRIX(:,9) * L * 1e-6];
        end
    end
    
    satisf_rand(k,:) = satisfaction_ratio{1}(k,:);
    satisf_rssi(k,:) = satisfaction_ratio{2}(k,:);
    satisf_lr(k,:) = satisfaction_ratio{3}(k,:);
    satisf_nn(k,:) = satisfaction_ratio{4}(k,:);
    
end

%% PLOTS
figure
subplot(2,1,1)
hAxes = gca;
x1 = satisf_rssi(:,1:8);
avg_x1 = mean(x1(:));
x2 = satisf_nn(:,1:8);
avg_x2 = mean(x2(:));
x3 = satisf_rssi(:,9:16);
avg_x3 = mean(x3(:));
x4 = satisf_nn(:,9:16);
avg_x4 = mean(x4(:));
x5 = satisf_rssi(:,17:24);
avg_x5 = mean(x5(:));
x6 = satisf_nn(:,17:24);
avg_x6 = mean(x6(:));
x = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:)];
boxplot(x,'colors',repmat('rb',1,3))
hold on
plot([2.5 2.5],[0 1.1],'k','linewidth',2)
plot([4.5 4.5],[0 1.1],'k','linewidth',2)
plot([avg_x1 avg_x2 avg_x3 avg_x4 avg_x5 avg_x6],'p','MarkerEdgeColor','k','MarkerFaceColor','yellow','MarkerSize',15);
title('Satisfaction (load vs throughput)')
ylabel('Satisfaction ratio')
xlabel('No. of STAs')
xticks([1.5 3.5 5.5])
xticklabels({'1-8', '9-16', '17-24'})
%legend('SSF', 'Vanila NN', 'Average');
set(gca,'fontsize',16)
grid on

subplot(2,1,2)
x1 = 1e-6.*min_tpt{1}(2,:)';
x2 = 1e-6.*min_tpt{1}(4,:)';
x3 = 1e-6.*min_tpt{2}(2,:)';
x4 = 1e-6.*min_tpt{2}(4,:)';
x5 = 1e-6.*min_tpt{3}(2,:)';
x6 = 1e-6.*min_tpt{3}(4,:)';
avg_x1 = mean(x1);
avg_x2 = mean(x2);
avg_x3 = mean(x3);
avg_x4 = mean(x4);
avg_x5 = mean(x5);
avg_x6 = mean(x6);
x = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:)];
boxplot(x,'colors',repmat('rb',1,3))
hold on
plot([2.5 2.5],[0 1.1*max(max(x))],'k','linewidth',2)
plot([4.5 4.5],[0 1.1*max(max(x))],'k','linewidth',2)
plot([avg_x1 avg_x2 avg_x3 avg_x4 avg_x5 avg_x6],'p','MarkerEdgeColor','k','MarkerFaceColor','yellow','MarkerSize',15);
title('Minimum throughput')
ylabel('Min. throughput [Mbps]')
xlabel('No. of STAs')
xticks([1.5 3.5 5.5])
xticklabels({'1-8', '9-16', '17-24'})
legend({'SSF', 'Vanilla NN', 'Average'},'NumColumns',3);
set(gca,'fontsize',16)
grid on

%%
figure
x1 = 1e-6.*tpt{1}(2,:,:);
x2 = 1e-6.*tpt{1}(4,:,:);
x3 = 1e-6.*tpt{2}(2,:,:);
x4 = 1e-6.*tpt{2}(4,:,:);
x5 = 1e-6.*tpt{3}(2,:,:);
x6 = 1e-6.*tpt{3}(4,:,:);
avg_x1 = mean(mean(x1));
avg_x2 = mean(mean(x2));
avg_x3 = mean(mean(x3));
avg_x4 = mean(mean(x4));
avg_x5 = mean(mean(x5));
avg_x6 = mean(mean(x6));
x = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:)];
boxplot(x,'colors',repmat('rb',1,3))
hold on
plot([2.5 2.5],[0 1.1*max(max(x))],'k','linewidth',2)
plot([4.5 4.5],[0 1.1*max(max(x))],'k','linewidth',2)
plot([avg_x1 avg_x2 avg_x3 avg_x4 avg_x5 avg_x6],'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
title('Throughput')
ylabel('Throughput [Mbps]')
xlabel('No. of STAs')
xticks([1.5 3.5 5.5])
xticklabels({'1-8', '9-16', '17-24'})
legend({'SSF', 'Vanilla NN', 'Average'},'NumColumns',3);
set(gca,'fontsize',16)
grid on