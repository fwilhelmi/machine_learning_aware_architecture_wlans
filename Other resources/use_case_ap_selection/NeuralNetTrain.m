%% Load Dataset
clear
close all
clc

filename = 'output_stas.csv';

% Read the .csv containing the dataset
T = readtable(filename);
T = T{:,:};
L = 12000;  % Fixed packet size [b]

% Convert table to arrays
% X = [T(:,2)/90 T(:,4)/144.4 T(:,5) T(:,6)/1000];
figure
sgtitle('Empirical PDF of the input features')

RSSI_dbm = T(:,2);  % [dBm]
RSSI_mW = 10.^(RSSI_dbm/10);
RSSI_mW_scaled = RSSI_mW / max(RSSI_mW);
subplot(2,3,1)
[f,xi] = ksdensity(RSSI_mW_scaled); 
plot(xi,f);
ylabel('P(RSSI)')
xlabel('RSSI')

rate_bps = T(:,3);  % [bps]
rate_bps_scaled = rate_bps / (144e6);
[f,xi] = ksdensity(rate_bps_scaled); 
subplot(2,3,2)
plot(xi,f);
ylabel('P(r)')
xlabel('r')

load_sta = T(:,4);  % [pkt/s]
load_sta_scaled = load_sta / max(load_sta);
[f,xi] = ksdensity(load_sta_scaled); 
subplot(2,3,3)
plot(xi,f);
ylabel('P(load_{STA})')
xlabel('load_{STA}')

deliv_ratio = T(:,5);   % unitless
deliv_ratio_scaled = deliv_ratio / max(deliv_ratio);
subplot(2,3,4)
[f,xi] = ksdensity(deliv_ratio_scaled); 
plot(xi,f);
ylabel('P(deliv\_ratio)')
xlabel('deliv\_ratio')

load_ap = T(:,6);   % [pkt/s]
load_ap_scaled = load_ap / max(load_ap);
subplot(2,3,5)
[f,xi] = ksdensity(load_ap_scaled); 
plot(xi,f);
ylabel('P(load_{AP})')
xlabel('load_{AP}')

% Output
throughput_sta = T(:,8);    % [pkt/s]
%throughput_sta_bps = throughput_sta * L * 1e-6;
throughput_sta_scaled = throughput_sta*L / 10*1e-6;
[f,xi] = ksdensity(throughput_sta_scaled); 
subplot(2,3,6)
plot(xi,f);
ylabel('P(thr)')
xlabel('thr')

%%  TRAIN ARTIFICIAL NEURAL NETWORK

fprintf("\n\n*** TRAIN ARTIFICIAL NEURAL NETWORK ***\n")

num_samples = length(rate_bps_scaled);

X = [rate_bps_scaled load_sta_scaled deliv_ratio_scaled load_ap_scaled]';
X_train = X(:,1:(num_samples/2));
X_test = X(:,(num_samples/2 + 1):num_samples);

Y = throughput_sta_scaled';
Y_train = Y(1:(num_samples/2));
Y_test = Y((num_samples/2 + 1):num_samples);

% no. of hidden layers
net = feedforwardnet([8 8]);
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
fprintf("- MAE: %f\n", mae)
fprintf("- MAE: %f Mbps\n", mae * max(throughput_sta) * L * 1e-6)

Y_pred_test = net(X_test);
mse = immse(Y_pred_test,Y_test);
mae = sqrt(mse);
fprintf("Test\n")
fprintf("- MSE: %f\n", mse)
fprintf("- MAE: %f\n", mae)
fprintf("- MAE: %f Mbps\n", mae * max(throughput_sta) * L * 1e-6)

% %%
% load('nn_function.mat')
% X = [0.1 0.2 0.3 0.4 0.5]';
% caro = net(X);
% 
%% 
save('nn_function.mat', 'net');