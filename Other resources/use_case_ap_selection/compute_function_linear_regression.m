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

filename = 'output_stas_2.csv';

% Read the .csv containing the dataset
T = readtable(filename);
T = T{:,:};
% Convert table to arrays
% X = [T(:,2)/90 T(:,4)/144.4 T(:,5) T(:,6)/1000];
X1 = [T(:,2)/90 T(:,3)/144.4 T(:,4)/833 T(:,5) T(:,6)/1000];
X2 = [T(:,3)/144.4 T(:,4)/833 T(:,5) T(:,6)/1000];
X3 = [T(:,2)/90 T(:,4)/833 T(:,5) T(:,6)/1000];
X4 = [T(:,2)/90 T(:,3)/144.4 T(:,5) T(:,6)/1000];
X5 = [T(:,2)/90 T(:,3)/144.4 T(:,4)/833 T(:,5)];
Y = T(:,8)/833;

% Compute linear regression
L1 = stepwisefit(X1,Y,'inmodel',true(1,5));
L2 = stepwisefit(X2,Y,'inmodel',true(1,4));
L3 = stepwisefit(X2,Y,'inmodel',true(1,4));)
L4 = stepwisefit(X4,Y,'inmodel',true(1,4));
L5 = stepwisefit(X5,Y,'inmodel',true(1,4));

% % PLOTS
% labels = {'RSSI' 'Load STA' 'Delivery ratio' 'Load AP'}; 
% data = [T(:,2)/90 T(:,4)/144.4 T(:,5) T(:,6)/1000];  
% [h,ax] = plotmatrix(data);                        
% % create a 4 x 4 matrix of plots 
% for i = 1:4                                       
%     % label the plots   
%     xlabel(ax(4,i), labels{i})   
%     ylabel(ax(i,1), labels{i}) 
% end