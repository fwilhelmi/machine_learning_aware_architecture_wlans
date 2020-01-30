%-------------------------------------------------------------------------
% WIFIX Platform - Channel load aware AP/Extender selection mechanism
%-------------------------------------------------------------------------
% PlotTopologyExtenders.m --> Visual display of the resulting topology
%                             of Extenders
%-------------------------------------------------------------------------

function PlotTopologyExtenders(M,Pt,Sens,margin_R,f_backbone,PL_backbone,T_avg_R)

%Sizes
sz = size(M);
dev = sz(1);

% Plots
figure

for i=1:dev
    if i==1
        labels = {'AP'};
        scatter(M(i,1),M(i,2),10,'filled');
        text(M(i,1),M(i,2),labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
        xlim([-150 150]);
        ylim([-150 150]);
        hold;
    else
        labels = i-1;
        scatter(M(i,1),M(i,2),10,'filled');
        text(M(i,1),M(i,2),['R',num2str(labels)],'VerticalAlignment','bottom','HorizontalAlignment','right');
    end
    
end

title('Extender Topology');

for i=1:dev
    index = M(i,5);
    
    if (index ~= 0)
        xx = [M(i,1) M(index,1)];
        yy = [M(i,2) M(index,2)];
        
        if (index == 1)
            line(xx,yy,'Color','black','LineStyle','--');
        else
            line(xx,yy,'Color','red','LineStyle','--');
        end
    end
end

Lmax = Pt - Sens;
Dmax = MaxDistance(Lmax,f_backbone,PL_backbone);
circle(M(1,1),M(1,2),Dmax);

txt = '\rightarrow AP Coverage';
text(Dmax, 0, txt);

Lmargin = Pt - margin_R;
Dmargin = MarginDistance(Lmargin,f_backbone,PL_backbone);
circle(M(1,1),M(1,2),Dmargin);

txt = '\rightarrow Direct Connection';
text(Dmargin, 0, txt);

txt = ['Average Delay = ', num2str(T_avg_R*1000),' ms'];
text(60, -90, txt);

end

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end
