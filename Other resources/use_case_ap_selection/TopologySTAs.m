%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% TopologySTAs.m --> Computation of the resulting topology of STAs
%-------------------------------------------------------------------------
%   1) This algorithm shall be performed by the RRM.
%   2) The algorithm shall be conducted with the information extracted
%       from beacons of the AP and the already existing Extenders.
%    	Namely, these are the required fields:
%       -- RSSI (AP->STA)
%       -- RSSI (Extenders->STA)
%       -- #hops to the AP from each AP/Extender
%       -- #children per device (*)
%       -- Occupied airtime (in %) of the AP/Extender in the channel 
%           corresponding to the access link
%       -- Occupied airtime (in %) of the AP/Extender in the channel 
%           corresponding to the backbone link   
%-------------------------------------------------------------------------
%   * The result is a coverage map that shows through what AP/Extender
%     an STA would connect to the network.
%-------------------------------------------------------------------------

function [M,N,routing_table,S_STA,D_STA,U_STA,A_STA,S_R,D_R,U_R,A_R, NEW_MATRIX] = TopologySTAs(net,...
        map_STA,M,N,WIFI_std,f_backbone,f_access,PL_access,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,score_mode,max_STA_per_R,w_a,w_b,w_c,channel_load_ext)

% (Note that the first row of M corresponds to the AP, filled with 0s)
% M = [posX    posY    #hops   #Children_R  Parent_Index    Delay   Rate   Channel     #Children_S   lambda_R   DBPS    Airtime access  Airtime backbone]
%      1       2       3       4            5               6       7      8           9             10         11      12              13

% N = [posX    posY    Parent_Index    Channel     Rate    Lambda  DBPS     Type    Score_mode]
%      1       2       3               4           5       6       7        8       9     

%Parent_Index = 0: no connection
%Parent_Index = 1: connected to the AP
%Parent_Index = 2+: connected to Extender k

dev = size(M,1);
sta = size(N,1);

Dist_STA = zeros(1,sta);
PL_STA = zeros(1,sta);
Pr_STA = zeros(1,sta);

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

NEW_MATRIX = zeros(sta, 7);

% Radio Parameters
NC = 50000;       %Score of non-decoded packets (i.e., below Sens. level)

% Mechanism parameters
if (map_STA >= 1)
    c_low = -100;
    c_high = 100;
    res_plot = 0.5;
    [X,Y] = meshgrid(c_low:res_plot:c_high);
    Dev_chosen = zeros(length(X),length(Y));
else
    RSSI_map = zeros(1,dev);
end

% One iteration per each new STA
for m=1:sta        
    if (map_STA >= 1)  %Map method              
        % Creating values of the map (IGNORE)
        for i=1:length(X)
            for j=1:length(Y)
                for k=1:dev
                    D = sqrt((X(i,j)-M(k,1)).^2+(Y(i,j)-M(k,2)).^2);
                    PL = PathLossModel(f_access,D,PL_access);
                    RSSI_map(k) = Pt - PL;
                end
                Dev_chosen(i,j) = APExtenderSelector(M,N,m,f_access,RSSI_map,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,score_mode,max_STA_per_R,w_a,w_b,w_c);
            end
        end
        
        % Ploting intermediate maps
        if (map_STA >= 2)
            PlotTopologySTAs(M,N,X,Y,Dev_chosen);
        end
                
        % Allocating a new STA
        N(m,3) = Dev_chosen(N(m,2)*(1/res_plot)+(c_high-c_low)+1,N(m,1)*(1/res_plot)+(c_high-c_low)+1);
    else    %Fastest computation method
        for k=1:dev
            D = sqrt((N(m,1)-M(k,1)).^2+(N(m,2)-M(k,2)).^2);
            PL = PathLossModel(f_access,D,PL_access);
            RSSI_map(k) = Pt - PL;
        end
        [N(m,3), predicted_tpt] = APExtenderSelector(net,...
            M,N,m,f_access,RSSI_map,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,N(m,9),max_STA_per_R,w_a,w_b,w_c);
    end
    
    if (N(m,3) ~= 0)                %The STA has found a parent (1:AP/Other:Extender)
        Dist_STA(m) = sqrt((N(m,1)-M(N(m,3),1)).^2+(N(m,2)-M(N(m,3),2)).^2);
        PL_STA(m) = PathLossModel(f_access,Dist_STA(m),PL_access);
        Pr_STA(m) = Pt - PL_STA(m);
        N(m,4) = M(N(m,3),8);           %Setting STA channel
        [N(m,7),N(m,5)] = RatesWIFI(Pr_STA(m),Sens,f_access);
        M(N(m,3),9) = M(N(m,3),9) + 1;  %Adding a child to the parent
        %%%%%%%% NEW - PRINT INFO TO OUTPUT TEXT (RSSI; LOAD AP; LOAD STA)
        %%%%%%%% - Then, add the throughput (to be computed)
        NEW_MATRIX(m,1) = N(m,3);    % AP #
        NEW_MATRIX(m,2) = Pr_STA(m);        % RSSI [dBm]
        NEW_MATRIX(m,3) = N(m,5);           % RATE_STA [bps]
        NEW_MATRIX(m,4) = N(m,6);           % LOAD_STA [pkt/s]
        active_STA_AP = [];
        active_STA_AP = find((N(:,3) == N(m,3)));
        if (M(N(m,3), 10) == 0)
            NEW_MATRIX(m,5) = 1;
        else
            NEW_MATRIX(m,5) = M(N(m,3), 10) / (sum(N(active_STA_AP,6))-N(m,6));    % RATIO_DELIVER_SUCCESFUL_AP [%]
        end
        NEW_MATRIX(m,6) = M(N(m,3), 10);    % LOAD AP[pkt/s]
        NEW_MATRIX(m,7) = M(N(m,3), 12);    % AT_AP [Ratio]
    end
    
%     disp('Parameters STA:')
%     disp(['   - Selected AP = ' num2str(NEW_MATRIX(m,1))]) 
%     disp(['   - RSSI = ' num2str(NEW_MATRIX(m,2))]) 
%     disp(['   - RATE_STA = ' num2str(NEW_MATRIX(m,3))]) 
%     disp(['   - LOAD_STA = ' num2str(NEW_MATRIX(m,4))]) 
%     disp(['   - RATIO_DELIVER_SUCCESFUL_AP = ' num2str(NEW_MATRIX(m,5))]) 
%     disp(['   - Load AP = ' num2str(NEW_MATRIX(m,6))])     
    
    % Computing number of different channels
    subset_ch = find(N(:,4) > 0);
    dif_ch = transpose(unique(N(subset_ch,4)));
    n_ch = length(dif_ch);
    A_CH = zeros(1,n_ch);
    
    % Computing WIFIModel per each channel [STA - Extender]
    for i=1:n_ch
        active_STA_CH = [];
        active_STA_CH = find((N(:,3) ~=0) & (N(:,4) == dif_ch(i)));
        if (isempty(active_STA_CH) == 0)
            [S_STA(active_STA_CH),D_STA(active_STA_CH),U_STA(active_STA_CH),A_STA(active_STA_CH)] = WIFIModel(N(active_STA_CH,6),L,N(active_STA_CH,7),N(active_STA_CH,5),f_access,WIFI_std);
            A_CH(i) = sum(A_STA(active_STA_CH)) + channel_load_ext(dif_ch(i));               
            for j=1:dev
                if (dif_ch(i) == M(j,8))
                    M(j,12) = A_CH(i);
                end
            end
        end
    end
               
    if S_STA == 0
        NEW_MATRIX(m,8) = 0;    % Throughput [bps]
    else 
        NEW_MATRIX(m,8) = S_STA(m)/L;    % Throughput [pkt/s]
    end
    
    if score_mode == 101
        disp(['Predicted value: ' num2str(predicted_tpt(N(m,3)-1))])
        disp(['Actual value: ' num2str((NEW_MATRIX(m,8))')])
    end
    
    NEW_MATRIX(m,9) = predicted_tpt(N(m,3)-1) - (NEW_MATRIX(m,8)*L*1e-6);
    
    %Computing resulting generation rate per each Extender in function of STAs
    if (N(m,3) > 0)
        M(N(m,3),10) = M(N(m,3),10) + S_STA(m)/L;
    end
    
    %Sum all generation rates from Extenders at ring #2 to Extenders at ring #1
    %TO DO: Improve this
    if (N(m,3) > 1)             % STA m is connected to a Extender
        if (M(N(m,3),3) == 2)   % This Extender is 2 hops away from the AP
            M(M(N(m,3),5),10) = M(M(N(m,3),5),10) + S_STA(m)/L;
        end
    end
    
    %WIFI Model for Extenders directly connected to the AP (1 hop)
    R_hop1 = find((M(:,3) == 1));
    if (isempty(R_hop1) == 0)               
        [S_R(R_hop1),D_R(R_hop1),U_R(R_hop1),A_R(R_hop1)] =  WIFIModel(M(R_hop1,10),L,M(R_hop1,11),M(R_hop1,7),f_backbone,WIFI_std);
        A_R_T = sum(A_R(R_hop1));        
        
        if (A_R_T > 1)
            A_R_T = 1;
        end
        
        M(R_hop1,13) = A_R_T;
    end
    
    %WIFI Model for Extenders at ring #2
    %TO DO: Improve this computation as it is not accurate enough
    for i=2:dev
        if (M(i,3) == 2)
            [S_R_aux,D_R_aux,U_R_aux,A_R_aux] = WIFIModel(M(i,10),L,M(i,11),M(i,7),f_backbone,WIFI_std);        
            M(i,13) = A_R_aux;
        end
    end
        
    %Updating satisfaction of Extenders
    R_loaded = find((M(2:dev,10) > 0)) + 1;
    SS_R(R_loaded) = S_R(R_loaded)./(transpose(M(R_loaded,10)).*L);
end

%Fill index_table
%The first row corresponds to the AP. The rest of rows are the Extenders
%The first column corresponds to the channel employed in the connection with STAs
routing_table(:,1) = M(:,8);
for i=1:sta
    if (N(i,3) ~= 0)
        vv = find(routing_table(N(i,3),:) == 0);
        routing_table(N(i,3),vv(1)) = i;
    end
end

%Set associated delay to Extenders at ring #2
%TO DO: Improve this computation as it is not accurate enough
for i=2:dev
    if (M(i,3) == 2)
        D_R(i) = M(i,6) + D_R(M(i,5));
    end
end

%Adding delay from Extenders to STAs
for i=1:sta
    if (N(i,3) ~= 0)
        D_STA(i) = D_STA(i) + D_R(N(i,3));
    end
end

% Update cumulated packets in the AP with traffic from Extenders
R_hop1 = transpose(find((M(:,3) == 1)));
S_R_AP = sum(S_R(R_hop1));
M(1,10) = M(1,10) + S_R_AP/L;

% Last iteration to print map
if (map_STA >= 1)
    for i=1:length(X)
        for j=1:length(Y)
            for k=1:dev
                D = sqrt((X(i,j)-M(k,1)).^2+(Y(i,j)-M(k,2)).^2);
                PL = PathLossModel(f_access,D,PL_access);
                RSSI_map(k) = Pt - PL;
            end
            Dev_chosen(i,j) = APExtenderSelector(M,f_access,RSSI_map,Pt,Sens,L,TPHY,SIFS,DIFS,Tslot,score_mode,max_STA_per_R,w_a,w_b,w_c);
        end
    end
       
    PlotTopologySTAs(M,N,X,Y,Dev_chosen);
end
end
