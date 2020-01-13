%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% PosGeneratorCircle.m --> This function places n_STA stations into an area 
%                    limited by the radius R and with the center 
%                    established in the position pos_AP
%-------------------------------------------------------------------------

function [pos_STA] = PosGeneratorCircle(center,n_STA,R,margin_over,margin_under)

pos_STA = zeros(n_STA,2);
a = margin_under/100;
b = 1+(margin_over/100);

    for i=1:n_STA
        seed1 = a + (b-a)*rand;
        seed2 = rand;

        pos_aux1 = seed1*R*cos(seed2*2*pi);
        pos_STA(i,1) = pos_aux1 + center(1);

        pos_aux2 = seed1*R*sin(seed2*2*pi);
        pos_STA(i,2) = pos_aux2 + center(2);
    end
end

%OLD VERSION:
%(To be used when 'map_STA = 1')
%
% function [pos_STA] = PosGeneratorCircle(center,n_STA,R,margin_over,margin_under)
% 
% pos_STA = zeros(n_STA,2);
% a = margin_under/100;
% b = 1+(margin_over/100);
% 
%     for i=1:n_STA
%         seed1 = a + (b-a)*rand;
%         seed2 = rand;
% 
%         pos_aux1 = seed1*R*cos(seed2*2*pi);
%         if (pos_aux1 > 0)
%             pos_STA(i,1) = floor(pos_aux1);
%         else
%             pos_STA(i,1) = ceil(pos_aux1);
%         end
%         pos_STA(i,1) = pos_STA(i,1) + center(1);
%         pos_STA(i,1) = pos_aux1 + center(1);
% 
%         pos_aux2 = seed1*R*sin(seed2*2*pi);
%         if (pos_aux2 > 0)
%             pos_STA(i,2) = floor(pos_aux2);
%         else
%             pos_STA(i,2) = ceil(pos_aux2);
%         end
%         pos_STA(i,2) = pos_STA(i,2) + center(2);
%         pos_STA(i,2) = pos_aux2 + center(2);
%     end
% end