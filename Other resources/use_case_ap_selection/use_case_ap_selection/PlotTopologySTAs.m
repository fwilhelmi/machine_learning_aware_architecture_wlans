%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% PlotTopologySTAs.m --> Visual display of the resulting topology of STAs
%-------------------------------------------------------------------------

function PlotTopologySTAs(M,N,X,Y,Dev_chosen)

% Parameters
dev = size(M,1);
sta = size(N,1);

% Plots
figure

%surf(X,Y,Z);
contourf(X,Y,Dev_chosen);
hold;
for i=1:dev
    if i==1
        labels2 = M(i,8);
        scatter3(M(i,1),M(i,2),10,'filled');
        text(M(i,1),M(i,2),['AP (CH #',num2str(labels2),')'],'VerticalAlignment','bottom','HorizontalAlignment','right');
    else
        labels = i-1;
        labels2 = M(i,8);
        scatter3(M(i,1),M(i,2),10,'filled');
        text(M(i,1),M(i,2),['R',num2str(labels),' (CH #',num2str(labels2),')'],'VerticalAlignment','bottom','HorizontalAlignment','right');
    end
end

for i=1:sta
    labels = i;
    scatter3(N(i,1),N(i,2),10,'filled');
    text(N(i,1),N(i,2),['STA',num2str(labels)],'VerticalAlignment','bottom','HorizontalAlignment','right');
end

for i=1:dev
    index = M(i,5);
    
    if (index > 0)
        xx = [M(i,1) M(index,1)];
        yy = [M(i,2) M(index,2)];
        
        if (index == 1)
            line(xx,yy,'Color','black','LineStyle','--');
        else
            line(xx,yy,'Color','red','LineStyle','--');
        end
    end
end

title('STA Topology');

for i=1:sta
    index = N(i,3);
    
    if (index > 0)
        xx = [N(i,1) M(index,1)];
        yy = [N(i,2) M(index,2)];
        line(xx,yy,'Color','blue','LineStyle','--');
    end
end

end