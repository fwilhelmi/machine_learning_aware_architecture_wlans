%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% MaxDistance.m --> Returns the maximum distance covered by a signal
%                   according to its physical parameters
%-------------------------------------------------------------------------

function Dmax = MaxDistance(PLmax,fc,PL_model)

switch PL_model
    
    case 1                              %ITU-R model
        Lo = 20*log10(fc/1E6)-28;           %Path loss intercept
        alfa = 3.1;                         %Attenuation factor
        
    case 2                              %TMB model (does not work!)
        Lo = 54.12;                         %Path loss intercept
        alfa = 2.06067;                     %Attenuation factor
        
    case 3                              %IEEE 802.11ax Residential
        Lo = 40.05 + 20*log10(fc/2.4E9) + 20*log10(5);
        alfa = 3.5;
        
    case 4                              %IEEE 802.11ax Enterprise
        Lo = 40.05 + 20*log10(fc/2.4E9) + 20*log10(10);
        alfa = 3.5;
        
    otherwise                           %Boris model (¿?)
        Lo = 30;                            %Path loss intercept
        alfa = 4.2;                         %Attenuation factor
end

Dmax = 10^((PLmax-Lo)/(10*alfa));

end