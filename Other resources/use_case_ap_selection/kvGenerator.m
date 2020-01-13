%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% kvGenerator.m --> This function gives kv_share (%) of STAs allocated 
%                   with the IEEE 802.11k/v capabilities
%-------------------------------------------------------------------------

function [type_STA,score_mode_STA] = kvGenerator(sta,kv_share,score_mode)

type_STA = zeros(sta,1);
score_mode_STA = zeros(sta,1);

for i = 1:sta
    
    if (rand()*100 < kv_share)
        type_STA(i) = 1;
        score_mode_STA(i) = score_mode;
    end
    
end

end