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

filename = 'output_stas.csv';

% Read the .csv containing the dataset
T = readtable(filename);
T = T{:,:};
% Convert table to arrays
X = [T(:,2)/90 T(:,4)/144.4 T(:,5) T(:,6)/1000];
Y = T(:,8)/833;

% Add padding (1s) to the dataset
%X = [ones(1, size(x,1))' x]; 
% Compute linear regression
%b = regress(y, X);
stepwisefit(X,Y,'inmodel',true(1,4))

% PLOTS
labels = {'RSSI' 'Load STA' 'Delivery ratio' 'Load AP'}; 
data = [T(:,2)/90 T(:,4)/144.4 T(:,5) T(:,6)/1000];  
[h,ax] = plotmatrix(data);                        
% create a 4 x 4 matrix of plots 
for i = 1:4                                       
    % label the plots   
    xlabel(ax(4,i), labels{i})   
    ylabel(ax(i,1), labels{i}) 
end