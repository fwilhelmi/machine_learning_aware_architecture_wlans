%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% MarginDistance.m --> Auxiliar function to plot resulting topology of
%                      Extenders                       
%-------------------------------------------------------------------------

function Dmargin = MarginDistance(Lmargin,fc,PL_model)

switch PL_model
    case 2                              %ITU-R model
        Lo = 20*log10(fc/1E6)-28;           %Path loss intercept
        alfa = 3.1;                         %Attenuation factor
    case 3                              %TMB model
        Lo = 54.12;                         %Path loss intercept
        alfa = 2.06067;                     %Attenuation factor
    otherwise                           %Boris model (¿?)
        Lo = 30;                            %Path loss intercept
        alfa = 4.2;                         %Attenuation factor
end

Dmargin = 10^((Lmargin-Lo)/(10*alfa));

end