function [ graf ] = geometry_cross( x_t,x_r,y_t,y_r,x_SD,y_SD,x_MD, y_MD,x_DI, y_DI, l_road, W_road, N_lanes, R_ce, R_ci, environment, scenario,Q,R )
%GEOMETRY Summary of this function goes here
%   Detailed explanation goes here
figure(1)
graf=x_t;

if environment==1 % Highway
    if scenario==0 % Highway Straight
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
    else % Highway crossing

        x=[0,l_road]; % lines
        plot (x,[0 0],'-k','LineWidth',2)
        hold on
        plot (x,[-W_road/N_lanes -W_road/N_lanes], '-k')
        hold on
        x1=[sqrt((R_ci.^2)-(R_ce.^2)),(l_road)-sqrt((R_ci.^2)-(R_ce.^2))]; % lines
        plot (x1,[-2*W_road/N_lanes -2*W_road/N_lanes], '-k')
        hold on
        plot (x,[W_road/N_lanes W_road/N_lanes], '-k')
        hold on
        plot (x1,[2*W_road/N_lanes 2*W_road/N_lanes], '-k')
        hold on
        % crossing down part
        theta_max2=asin(R_ce/R_ci)*360/(2*pi); % DEG
        theta_l1=0:0.01:pi/2;
        theta_l2=0:0.01:theta_max2*2*pi/360;
        
        [X_l1,Y_l1]=pol2cart(theta_l2,R_ci);
        plot (X_l1,-(W_road/2+R_ce)+Y_l1, '-k')
        hold on
        [X_l2,Y_l2]=pol2cart(theta_l1,R_ce);
        plot (X_l2,-(W_road/2+R_ce)+Y_l2, '-k')
        hold on
        [X_l3,Y_l3]=pol2cart(theta_l2+pi-theta_max2*2*pi/360,R_ci);
        plot (l_road+X_l3,-(W_road/2+R_ce)+Y_l3, '-k')
        hold on
        [X_l4,Y_l4]=pol2cart(theta_l1+pi/2,R_ce);
        plot (l_road+X_l4,-(W_road/2+R_ce)+Y_l4, '-k')
        hold on
                 
        % crossing up part
        [X_l1u,Y_l1u]=pol2cart(theta_l2+pi,R_ci);
        plot (l_road+X_l1u,(W_road/2+R_ce)+Y_l1u, '-k')
        hold on
        [X_l2u,Y_l2u]=pol2cart(theta_l1+pi,R_ce);
        plot (l_road+X_l2u,(W_road/2+R_ce)+Y_l2u, '-k')
        hold on
        [X_l3u,Y_l3u]=pol2cart(theta_l2+2*pi-theta_max2*2*pi/360,R_ci);
        plot (X_l3u,(W_road/2+R_ce)+Y_l3u, '-k')
        hold on
        [X_l4u,Y_l4u]=pol2cart(theta_l1+pi*3/2,R_ce);
        plot (X_l4u,(W_road/2+R_ce)+Y_l4u, '-k')
        hold on
    end
else % Rural
    if scenario==0 % Rural Straight
        plot ([0 l_road],[0 0],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([0 l_road],[4 4],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([0 l_road],[-4 -4],'LineWidth',2,'Color',[0 0 0]);
        hold on
    else % Rural crossing
        plot ([0 (l_road/2-W_road/2)],[0 0],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([(l_road/2+W_road/2) l_road],[0 0],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([0 (l_road/2-W_road/2)],[4 4],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([(l_road/2+W_road/2) l_road],[4 4],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([0 (l_road/2-W_road/2)],[-4 -4],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([(l_road/2+W_road/2) l_road],[-4 -4],'LineWidth',2,'Color',[0 0 0]);
        hold on
        
        plot ([l_road/2 l_road/2],[4 (l_road/2)],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([l_road/2 l_road/2],[(-l_road/2) -4],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([l_road/2-W_road/2 l_road/2-W_road/2],[4 (l_road/2)],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([l_road/2-W_road/2 l_road/2-W_road/2],[(-l_road/2) -4],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([l_road/2+W_road/2 l_road/2+W_road/2],[4 (l_road/2)],'LineWidth',2,'Color',[0 0 0]);
        hold on
        plot ([l_road/2+W_road/2 l_road/2+W_road/2],[(-l_road/2) -4],'LineWidth',2,'Color',[0 0 0]);
        hold on
    end
end

plot (x_t, y_t,'-bs','MarkerSize',20);
hold on
plot (x_r, y_r,'-rs','MarkerSize',20);


plot (x_SD(1:Q), y_SD,'md','MarkerSize',10);

hold on

plot (x_MD, y_MD,'ms','MarkerSize',20);

hold on

plot (x_DI(1:R), y_DI,'k.','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');

hold on

% theta_pat=0:0.1:2*pi;
% for i=1:length(x_t)
% [X_apt,Y_apt]=pol2cart(cos(theta_pat-phi_t-1.3),g_T);
% %x_tp=repmat(x_t,1,length(X_apt));
% %y_tp=repmat(y_t,1,length(Y_apt));
%
% plot (x_t(i)+X_apt,y_t(i)+Y_apt, '-g')
% hold on
%
% [X_apr,Y_apr]=pol2cart(cos(theta_pat-phi_r-1.3),g_R);
% % x_rp=repmat(x_r,1,length(theta_pat));
% % y_rp=repmat(y_r,1,length(theta_pat));
% plot (x_r(i)+X_apr,y_r(i)+Y_apr, '-g')
% hold on
% end
grid on
xlabel ('road lenght [m]');
ylabel ('road width [m]');
end

