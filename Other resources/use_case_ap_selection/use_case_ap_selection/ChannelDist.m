%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% ChannelDist.m --> Randomly distributes the available channels among the 
%                   existing devices
%-------------------------------------------------------------------------

function channel_R = ChannelDist(fc,dev)

if (fc == 2.4E9)
    
    channel_R = randi([1 3],1,dev);
    
    for i = 1:dev
        if (channel_R(i) == 2)
            channel_R(i) = 6;
        end
        if (channel_R(i) == 3)
            channel_R(i) = 11;
        end
    end
else
    disp('Error in ChannelDist function');
end

end