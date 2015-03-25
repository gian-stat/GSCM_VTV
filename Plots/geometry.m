function [ graf ] = geometry( x_t,x_r,y_t,y_r,x_SD,y_SD,x_MD, y_MD,x_DI, y_DI,Q,R, l_road, W_road, environment )
%GEOMETRY Summary of this function goes here
%   Detailed explanation goes here
figure(1)
graf=x_t;

if environment==1
    N_lanes=4;
    plot ([0 l_road],[0 0],'LineWidth',2,'Color',[0 0 0]);
    hold on
    plot ([0 l_road],[W_road/N_lanes W_road/N_lanes],'LineWidth',2,'Color',[0 0 0]);
    hold on
    plot ([0 l_road],[-W_road/N_lanes -W_road/N_lanes],'LineWidth',2,'Color',[0 0 0]);
    hold on
    plot ([0 l_road],[2*W_road/N_lanes 2*W_road/N_lanes],'LineWidth',2,'Color',[0 0 0]);
    hold on
    plot ([0 l_road],[-2*W_road/N_lanes -2*W_road/N_lanes],'LineWidth',2,'Color',[0 0 0]);
    hold on
else
    plot ([0 l_road],[0 0],'LineWidth',2,'Color',[0 0 0]);
    hold on
    plot ([0 l_road],[4 4],'LineWidth',2,'Color',[0 0 0]);
    hold on
    plot ([0 l_road],[-4 -4],'LineWidth',2,'Color',[0 0 0]);
    hold on
end

plot (x_t, y_t,'-bs','MarkerSize',20);
hold on
plot (x_r, y_r,'-rs','MarkerSize',20);

for j=1:ceil(Q)
plot (x_SD(j), y_SD(j),'md','MarkerSize',10);

hold on
end

plot (x_MD, y_MD,'ms','MarkerSize',20);

hold on

for j=1:ceil(R)
plot (x_DI(j), y_DI(j),'k.','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');

hold on
end

grid on
xlabel ('road lenght [m]');
ylabel ('road width [m]');
end

