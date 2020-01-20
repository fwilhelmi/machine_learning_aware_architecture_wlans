%% Load Dataset
clear
close all
clc

filename = 'output_stas.csv';

% Read the .csv containing the dataset
T = readtable(filename);
T = T{:,:};
num_original_samples = size(T,1);

%% DATASET PRUNING (NOT USED)
% stas_to_start = 6;
% stas_to_finish = 20;
% stas_comb_per_deploy = stas_to_finish - stas_to_start + 1;
% ratio = stas_comb_per_deploy / stas_to_finish;
% num_new_samples = ratio * num_original_samples;
% new_T = zeros(num_new_samples, 8);
% last_sample_ix_new = 0;
% fprintf('Pruning the dataset to have samples from %d to %d STAs...\n', stas_to_start, stas_to_finish)
% fprintf('- 1\n')
% for t = 1:20:num_original_samples
%     % fprintf('- %d\n', t)
%     if mod(t,10001) == 0
%         fprintf('- %.2f %%\n', t/num_original_samples * 100)
%     end   
%     first_sample_ix = t + stas_to_start - 1;
%     last_sample_ix = first_sample_ix - stas_to_start + stas_to_finish;
%     %fprintf('first_sample_ix = %d, last_sample_ix = %d\n', first_sample_ix, last_sample_ix)
%     info_to_add = T(first_sample_ix:last_sample_ix,:);
%     
%     first_sample_ix_new = last_sample_ix_new+1;
%     last_sample_ix_new = first_sample_ix_new + stas_to_finish - stas_to_start;
%     %fprintf('first_sample_ix_new = %d, last_sample_ix_new = %d\n', first_sample_ix_new, last_sample_ix_new)
%     new_T(first_sample_ix_new:last_sample_ix_new, :) = info_to_add;
% end
% 
% fprintf('Dataset pruned\n')

new_T = T;

%%

L = 12000;  % Fixed packet size [b]

% Convert table to arrays
% X = [T(:,2)/90 T(:,4)/144.4 T(:,5) T(:,6)/1000];
figure
sgtitle('Empirical PDF of the input features')

MIN_RSSI_DBM = -90;
RSSI_dbm = new_T(:,2);  % [dBm]
%RSSI_mW = 10.^(RSSI_dbm/10);
RSSI_scaled = 1 - RSSI_dbm / MIN_RSSI_DBM;
%RSSI_mW_scaled = RSSI_mW / max(RSSI_mW);
subplot(2,3,1)
[f,xi] = ksdensity(RSSI_scaled);
plot(xi,f);
ylabel('P(RSSI)')
xlabel('RSSI')

RATE_BPS_MAX = 144e6;
rate_bps = new_T(:,3);  % [bps]
rate_bps_scaled = rate_bps / RATE_BPS_MAX;
[f,xi] = ksdensity(rate_bps_scaled);
subplot(2,3,2)
plot(xi,f);
ylabel('P(r)')
xlabel('r')

LOAD_STA_MAX = 833;  % [pkt/s]
load_sta = new_T(:,4);  % [pkt/s]
load_sta_scaled = load_sta / LOAD_STA_MAX;
[f,xi] = ksdensity(load_sta_scaled);
subplot(2,3,3)
plot(xi,f);
ylabel('P(load_{STA})')
xlabel('load_{STA}')

DELIV_RATIO_MAX = 1;
deliv_ratio = new_T(:,5);   % unitless
deliv_ratio_scaled = deliv_ratio / DELIV_RATIO_MAX;
subplot(2,3,4)
[f,xi] = ksdensity(deliv_ratio_scaled);
plot(xi,f);
ylabel('P(deliv\_ratio)')
xlabel('deliv\_ratio')

LOAD_AP_MAX = 10 * LOAD_STA_MAX;
load_ap = new_T(:,6);   % [pkt/s]
load_ap_scaled = load_ap / LOAD_AP_MAX;
subplot(2,3,5)
[f,xi] = ksdensity(load_ap_scaled);
plot(xi,f);
ylabel('P(load_{AP})')
xlabel('load_{AP}')

% Output
THROUGHPUT_MAX = LOAD_STA_MAX;
throughput_sta = new_T(:,8);    % [pkt/s]
%throughput_sta_bps = throughput_sta * L * 1e-6;
throughput_sta_scaled = throughput_sta / THROUGHPUT_MAX;
[f,xi] = ksdensity(throughput_sta_scaled);
subplot(2,3,6)
plot(xi,f);
ylabel('P(thr)')
xlabel('thr')

%%  TRAIN ARTIFICIAL NEURAL NETWORK

fprintf("\n\n*** TRAIN ARTIFICIAL NEURAL NETWORK ***\n")

num_samples = length(rate_bps_scaled);

% Dataset
X = [rate_bps_scaled load_sta_scaled deliv_ratio_scaled load_ap_scaled]';
Y = throughput_sta_scaled';

% - Training dataset
X_train = X(:,1:(num_samples/2));
Y_train = Y(1:(num_samples/2));

% - Test dataset
X_test = X(:,(num_samples/2 + 1):num_samples);
Y_test = Y((num_samples/2 + 1):num_samples);

fprintf('No. of samples training dataset: %d\n',length(X_train))
fprintf('No. of samples testing dataset: %d\n',length(X_test))

% no. of hidden layers
net = feedforwardnet([16 8 4]);
%net.trainParam.mu_max = 1e10;
net.divideParam.trainRatio = 0.80; % training set [%]
net.divideParam.valRatio = 0.10; % validation set [%]
net.divideParam.testRatio = 0.10; % test set [%]
% train a neural network
net.trainParam.epochs = 30;
fprintf("Training ANN...\n")
[net,tr] = train(net,X_train,Y_train);
fprintf("- Training completed!\n")
%%
%---------------------------------
% view net
%view (net)
% simulate a network over complete input range
Y_pred_train = net(X_train);
mse = immse(Y_pred_train,Y_train);
mae = sqrt(mse);
fprintf("Train\n")
fprintf("- MSE: %f\n", mse)
fprintf("- MAE: %f [%%]\n", mae)
fprintf("- MAE: %f Mbps\n", mae * max(throughput_sta) * L * 1e-6)

Y_pred_test = net(X_test);
mse = immse(Y_pred_test,Y_test);
mae = sqrt(mse);
fprintf("Test\n")
fprintf("- MSE: %f\n", mse)
fprintf("- MAE: %f [%%]\n", mae)
fprintf("- MAE: %f Mbps\n", mae * max(throughput_sta) * L * 1e-6)

% %%
% load('nn_function.mat')
% X = [0.1 0.2 0.3 0.4 0.5]';
% caro = net(X);
%
%%
save('nn_function.mat', 'net');

%%

num_deploys_to_plot = 20;
figure

for ix_ventena = 1:num_deploys_to_plot
    subplot(num_deploys_to_plot,1,ix_ventena)
    mae = zeros(stas_comb_per_deploy,1);
    
    first_ix_sample = (ix_ventena-1) * stas_comb_per_deploy + 1;
    last_ix_sample = (ix_ventena-1) * stas_comb_per_deploy + stas_comb_per_deploy;
    
%     fprintf('first_ix_sample = %d, last_ix_sample = %d\n', first_ix_sample, last_ix_sample)
    
    for ix_sample = first_ix_sample : last_ix_sample
                
        x = X_train(:,ix_sample);
        y = Y_train(ix_sample);
        
        y_pred = net(x);
        mae(ix_sample - first_ix_sample + 1) = abs(y - y_pred);
%         fprintf("- Output (real/pred): %.2f /%.2f \n", y, y_pred)
%         fprintf("- Output (real/pred): %.2f /%.2f pkt/s\n", y * THROUGHPUT_MAX, y_pred *THROUGHPUT_MAX)
%         fprintf("- Output (real/pred): %.2f /%.2f Mbps\n", y * max(throughput_sta) * L * 1e-6, y_pred * max(throughput_sta) * L * 1e-6)
%         fprintf("- MAE: %f [%%]\n", mae(ix_sample))
%         fprintf("- MAE: %f Mbps\n", mae(ix_sample) * max(throughput_sta) * L * 1e-6)
        
    end
    plot(mae .* max(throughput_sta) * L * 1e-6)
    xticks({})
    ylabel(['# ' num2str(ix_ventena)])
end
xticks(1:15)
xticklabels(stas_to_start:stas_to_finish)
xlim([1 stas_comb_per_deploy])
sgtitle('MAE [Mbps]')
xlabel('# Last STA')

