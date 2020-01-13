%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% PosGeneratorRectangle.m --> This function places n_STA stations into a 
%                             rectangle defined by its center and the 
%                             length of its sides
%-------------------------------------------------------------------------

function [pos_STA] = PosGeneratorRectangle(center,n_STA,side_h,side_v)

pos_STA = zeros(n_STA,2);

    for i=1:n_STA
        
        pos_STA(i,1) = side_h*rand + center(1) - side_h/2;
        pos_STA(i,2) = side_v*rand + center(2) - side_v/2;
        
    end
end
