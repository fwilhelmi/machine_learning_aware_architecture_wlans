%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% PathLossModel.m --> Set of available path loss models together with 
%                     their physical equations according to the frequency
%                     fc and the distance d
%-------------------------------------------------------------------------

function PL = PathLossModel(fc,d,PL_model)

switch PL_model
    
    case 1                              %ITU-R model
        Lo = 20*log10(fc/1E6)-28;               %Path loss intercept
        alfa = 3.1;                             %Attenuation factor
        PL = Lo + 10.*alfa.*log10(d);
        
    case 2                              %TMB model (Only for 5 GHz)
        Lo = 54.12;                             %Path loss intercept
        alfa = 2.06067;                         %Attenuation factor
        k = 5.25;                               %Attenuation per wall
        W = 0.1467;                             %Averaged wall factor
        PL = Lo + 10.*alfa.*log10(d) + k*W*d;
        
    case 3                              %IEEE 802.11ax Residential
        if (d <= 5)
            Lo = 40.05 + 20*log10(fc/2.4E9);
            alfa = 5;
            k = 5;
            W = 0.1467;
            PL = Lo + 10.*alfa.*log10(d) + k*W*d;
        else
            Lo = 40.05 + 20*log10(fc/2.4E9) + 20*log10(5);
            alfa = 6.5;
            k = 5;
            W = 0.1467;
            PL = Lo + 10.*alfa.*log10(d/5) + k*W*d;
        end
              
%         n_walls = 1/10;
%         n_floors = 1/3;
%         L_iw = 5;
%         min_distance = d;
%         if d > 5, min_distance = 5; end
%         Lo = 40.05 + 20*log10(fc/2.4E9) + 20*log10(min_distance)...
%             + 18.3*((d*n_floors)^(((d*n_floors)+2)/((d*n_floors)+1)) ...
%             - 0.46) + L_iw*(d*n_walls);
%         if (d >= 5)
%             PL = Lo + 35*log10(d/5);
%         else
%             PL = Lo;
%         end
        
    case 4                              %IEEE 802.11ax Enterprise
        if (d <= 10)
            Lo = 40.05 + 20*log10(fc/2.4E9);
            alfa = 2;
            k = 7;
            W = 0.1467;
            PL = Lo + 10.*alfa.*log10(d) + k*W*d;
        else
            Lo = 40.05 + 20*log10(fc/2.4E9) + 20*log10(10);
            alfa = 3.5;
            k = 7;
            W = 0.1467;
            PL = Lo + 10.*alfa.*log10(d/10) + k*W*d;
        end
        
    otherwise                           %Boris model (�?)
        Lo = 30;                                %Path loss intercept
        alfa = 4.2;                             %Attenuation factor
        PL = Lo + 10.*alfa.*log10(d);
    
end

end