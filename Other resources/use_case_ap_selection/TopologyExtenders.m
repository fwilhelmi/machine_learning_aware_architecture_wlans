%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% TopologyExtenders.m --> Computation of the resulting topology of Extenders
%-------------------------------------------------------------------------
%   1) This algorithm shall be performed by the RRM.
%   2) The algorithm shall be conducted with the information extracted
%       from Performance Metrics Reports sent by the AP and Extenders.
%       Namely, these are the required fields:
%       -- RSSI (AP->newExtender)
%       -- RSSI (Extenders->newExtender) (all combinations)
%
%       RSSI is basically used to infer the used MCS in each combination
%       of a pair of devices and then to extract the time delay of a
%       packet transmission. For this reason, the algorithm could also be
%       run by the RRM if knowing any of these two metrics.
%   3) Algorithm is limited to 2-hop routes
%-------------------------------------------------------------------------
%   * The result is the routing topology between the AP and the Extenders
%-------------------------------------------------------------------------

function [M] = TopologyExtenders(map_R,M,f_backbone,PL_backbone,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,ext_conn_alg,margin_R,max_R_per_R)

% (Note that the first row of M corresponds to the AP, filled with 0s)
% M = [posX    posY    #hops   #Children_R  Parent_Index    Delay   Rate   Channel     #Children_S   lambda_R   DBPS    Airtime access  Airtime backbone]
%      1       2       3       4            5               6       7      8           9             10         11      12              13

% Matrix dimensions
dev = size(M,1);
DBPS_dir = zeros(1,dev-1);
rate_dir = zeros(1,dev-1);
P = zeros(1,dev-1);

% Obtaining Pr and Delay per each Extender when directly connected to the AP
for i=2:dev
    D = sqrt((M(1,1)-M(i,1)).^2+(M(1,2)-M(i,2)).^2);
    PL = PathLossModel(f_backbone,D,PL_backbone);
    Pr = Pt - PL;
    P(i-1) = Pr;
    [DBPS_dir(i-1),rate_dir(i-1)] = RatesWIFI(Pr,Sens,f_backbone);
    T = PacketDelay(L,TPHY,DBPS_dir(i-1),SIFS,DIFS,Tslot);
    M(i,6) = T;
    
    if ((Pr >= margin_R) && (M(1,4) < max_R_per_R))
        M(i,3) = 1;
        M(1,4) = M(1,4) + 1;                %Update AP children
        M(i,5) = 1;
        M(i,7) = rate_dir(i-1);
        M(i,11) = DBPS_dir(i-1);
    end
end

%Setting best parent for each Extender
for i=2:dev
    if (M(i,5) == 0)
        rate_alt = 0;
        DBPS_alt = 0;
        best_Pr = P(i-1);
        best_T = M(i,6);
        
        for j=2:dev
            if ((M(j,3) == 1)&&(j ~= i))
                D_new = sqrt((M(i,1)-M(j,1)).^2+(M(i,2)-M(j,2)).^2);
                PL_new = PathLossModel(f_backbone,D_new,PL_backbone);
                Pr_new = Pt - PL_new;
                [DBPS_new,rate_new] = RatesWIFI(Pr_new,Sens,f_backbone);
                T_aux = PacketDelay(L,TPHY,DBPS_new,SIFS,DIFS,Tslot);
                T_new = T_aux + M(j,6);
                
                if (ext_conn_alg == 0)  %Based on RSSI (FON method)
                    if (Pr_new > best_Pr)
                        pos_alt = j;
                        best_Pr = Pr_new;
                        best_T = T_new;
                        rate_alt = rate_new;
                        DBPS_alt = DBPS_new;
                    end
                else                    %Based on delay (UPF method)
                    if (T_new < best_T)
                        pos_alt = j;
                        best_T = T_new;
                        rate_alt = rate_new;
                        DBPS_alt = DBPS_new;
                    end
                end
            end
        end
        
        if (ext_conn_alg == 0)          %Based on RSSI (FON method)
            if ((best_Pr > P(i-1)) && (M(pos_alt,4) < max_R_per_R))
                M(i,3) = 2;
                M(pos_alt,4) = M(pos_alt,4) + 1;
                M(i,5) = pos_alt;
                M(i,6) = best_T;
                M(i,7) = rate_alt;
                M(i,11) = DBPS_alt;
            else
                if ((P(i-1) > Sens) && (M(1,4) < max_R_per_R))
                    M(i,3) = 1;
                    M(1,4) = M(1,4) + 1;    %Update AP children
                    M(i,5) = 1;
                    M(i,7) = rate_dir(i-1);
                    M(i,11) = DBPS_dir(i-1);
                end
            end
        else                            %Based on delay (UPF method)
            if ((best_T < M(i,6)) && (M(pos_alt,4) < max_R_per_R))
                M(i,3) = 2;
                M(pos_alt,4) = M(pos_alt,4) + 1;
                M(i,5) = pos_alt;
                M(i,6) = best_T;
                M(i,7) = rate_alt;
                M(i,11) = DBPS_alt;
            else
                if ((M(i,6) ~= Inf) && (M(1,4) < max_R_per_R))
                    M(i,3) = 1;
                    M(1,4) = M(1,4) + 1;    %Update AP children
                    M(i,5) = 1;
                    M(i,7) = rate_dir(i-1);
                    M(i,11) = DBPS_dir(i-1);
                end
            end
        end
    end
end

%Obtaining average Delay to the AP
vec = M(:,6);
T_avg_R = mean(vec(find((vec>0)&(vec<Inf))));

%Plot
if (map_R == 1)
    PlotTopologyExtenders(M,Pt,Sens,margin_R,f_backbone,PL_backbone,T_avg_R);
end

end
