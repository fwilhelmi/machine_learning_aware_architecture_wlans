%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% MeanGenerator.m --> Auxiliar function to create averages of network
%                     performance metrics               
%-------------------------------------------------------------------------

function [S_T_avg,E_T_avg,share_ok,D_avg_avg,D_max_avg,SS_avg_avg,SS_min_avg,assoc_STA_avg,assoc_STA_AP_avg,assoc_STA_E_avg] = MeanGenerator(S_T,E_T,E_T_ok,D_avg,D_max,SS_avg,SS_min,assoc_STA,assoc_STA_AP,assoc_STA_E,num_rep)

S_T_avg = mean(S_T);
E_T_avg = mean(E_T);
share_ok = (sum(E_T_ok)/num_rep)*100;
D_avg_avg = mean(D_avg);
D_max_avg = mean(D_max);
SS_avg_avg = mean(SS_avg);
SS_min_avg = mean(SS_min);
assoc_STA_avg = mean(assoc_STA);
assoc_STA_AP_avg = mean(assoc_STA_AP);
assoc_STA_E_avg = mean(assoc_STA_E);

end