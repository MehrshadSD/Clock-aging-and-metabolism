% script to set up common properties for panels

set(gca,'Xtick',0:6:48,'Fontsize',16)
set(gca,'FontName','Arial')
xlabel('Time (h)')
% semi-transparent rectangles during dark cycles
rectangle('Position',[12 0 12 100],'FaceColor',[.5 .5 .5 .2])
rectangle('Position',[36 0 12 100],'FaceColor',[.5 .5 .5 .2])
% set line appearance
c = linspecer(4);
l1.Color = c(1,:);
l1.LineWidth = 4;
l1.LineStyle = '-';
if exist('l2','var')
    l2.Color = c(2,:);
    l2.LineWidth = 4;
    l2.LineStyle = '-';
end
if exist('l3','var')
    l3.Color = c(3,:);
    l3.LineWidth = 4;
    l3.LineStyle = '-';
end
if exist('l4','var')
    l4.Color = c(4,:);
    l4.LineWidth = 4;
    l4.LineStyle = '-';
end
