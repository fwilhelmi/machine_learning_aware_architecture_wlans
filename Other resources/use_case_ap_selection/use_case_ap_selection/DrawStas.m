%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% PathLossModel.m --> Set of available path loss models together with 
%                     their physical equations according to the frequency
%                     fc and the distance d
%-------------------------------------------------------------------------

function [] = DrawStas(pos_STA, pos_R)

    plot(pos_STA(:,1), pos_STA(:,2), 'x', 'markersize', 10);
    hold on
    plot(pos_R(:,1), pos_R(:,2), 'o', 'markersize', 15);
    legend({'STA', 'AP'})
end