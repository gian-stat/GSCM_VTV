%% Environment 1: Rural straight
% You can change your table to structure like this
% Advantage: all values are associated with a name. Good readability
% But you may have to type all the values and change your code
% If you think, this is just too much, ignore and do not spend time on this

paratable(1).chi.MD = 0.001;
paratable(1).chi.SD = 0.05;
paratable(1).chi.DI = 1;

paratable(1).n.LOS  = 1.6;
paratable(1).n.MD   = 3.5*rand(1,ceil(paratable(1).chi.MD*l_road));
paratable(1).n.SD   = 3.5*rand(1,ceil(paratable(1).chi.SD*l_road));
paratable(1).n.DI   = 3.0;

paratable(1).G0.LOS = -9;
paratable(1).G0.MD  = -89+24*paratable(1).n.MD;
paratable(1).G0.SD  = -89+24*paratable(1).n.SD;
paratable(1).G0.DI  = -23;


paratable(1).mu_sigma.LOS = 11.7;
paratable(1).mu_sigma.MD  = 15.1;
paratable(1).mu_sigma.SD  = 14.8;

paratable(1).mu_c.LOS = 8;
paratable(1).mu_c.MD  = 8.3;
paratable(1).mu_c.SD  = 2.5;

paratable(1).d_c_min.LOS = 5.4;
paratable(1).d_c_min.MD  = 2.5;
paratable(1).d_c_min.SD  = 1.4;


paratable(1).y1.SD = -5.5;
paratable(1).y1.DI = -9.5;

paratable(1).y2.SD = 5.5;
paratable(1).y2.DI = 9.5;

paratable(1).W_DI    = 5;
paratable(1).W_road  = 8;
paratable(1).N_lanes = 2;

%% Environment 2: Highway straight

paratable(2).chi.MD = 0.005;
paratable(2).chi.SD = 0.005;
paratable(2).chi.DI = 1;

paratable(2).n.LOS  = 1.8;
paratable(2).n.MD   = 3.5*rand(1,ceil(paratable(2).chi.MD*l_road));
paratable(2).n.SD   = 3.5*rand(1,ceil(paratable(2).chi.SD*l_road));
paratable(2).n.DI   = 5.4;

paratable(2).G0.LOS = -5;
paratable(2).G0.MD  = -89+24*paratable(2).n.MD;
paratable(2).G0.SD  = -89+24*paratable(2).n.SD;
paratable(2).G0.DI  = -104;

paratable(2).mu_sigma.LOS = 6.8;
paratable(2).mu_sigma.MD  = 9.4;
paratable(2).mu_sigma.SD  = 6.3;

paratable(2).mu_c.LOS = 7.2;
paratable(2).mu_c.MD  = 5.4;
paratable(2).mu_c.SD  = 4.9;

paratable(2).d_c_min.LOS = 4.4;
paratable(2).d_c_min.MD  = 1.1;
paratable(2).d_c_min.SD  = 1.0;


paratable(2).y1.SD = -10;
paratable(2).y1.DI = -13.5;

paratable(2).y2.SD = 10;
paratable(2).y2.DI = 13.5;

paratable(2).W_DI    = 5;
paratable(2).W_road  = 18;
paratable(2).N_lanes = 4;

%% Environment 1: Rural intersection

paratable(3).chi.MD = 0.001;
paratable(3).chi.SD = 0.05;
paratable(3).chi.DI = 1;

paratable(3).n.LOS  = 1.6;
paratable(3).n.MD   = 3.5*rand(1,ceil(paratable(3).chi.MD*(l_road-8)*2));
paratable(3).n.SD   = 3.5*rand(1,ceil(paratable(3).chi.SD*(l_road-8)*2));
paratable(3).n.DI   = 3.0;

paratable(3).G0.LOS = -9;
paratable(3).G0.MD  = -89+24*paratable(3).n.MD;
paratable(3).G0.SD  = -89+24*paratable(3).n.SD;
paratable(3).G0.DI  = -23;

paratable(3).mu_sigma.LOS = 11.7;
paratable(3).mu_sigma.MD  = 15.1;
paratable(3).mu_sigma.SD  = 14.8;

paratable(3).mu_c.LOS = 8;
paratable(3).mu_c.MD  = 8.3;
paratable(3).mu_c.SD  = 2.5;

paratable(3).d_c_min.LOS = 5.4;
paratable(3).d_c_min.MD  = 2.5;
paratable(3).d_c_min.SD  = 1.4;

paratable(3).y1.SD = -5.5;
paratable(3).y1.DI = -9.5;

paratable(3).y2.SD = 5.5;
paratable(3).y2.DI = 9.5;

paratable(3).W_DI    = 5;
paratable(3).W_road  = 8;
paratable(3).N_lanes = 2;

%% Environment 2: Highway intersection

paratable(4).chi.MD = 0.05;
paratable(4).chi.SD = 0.05;
paratable(4).chi.DI = 1.5;

paratable(4).n.LOS  = 1.8;
paratable(4).n.MD   = 3.5*rand(1,ceil(paratable(4).chi.MD*l_road));
paratable(4).n.SD   = 3.5*rand(1,ceil(paratable(4).chi.SD*l_road));
paratable(4).n.DI   = 5.4;

paratable(4).G0.LOS = -5;
paratable(4).G0.MD  = -89+24*paratable(4).n.MD;
paratable(4).G0.SD  = -89+24*paratable(4).n.SD;
paratable(4).G0.DI  = -104;

paratable(4).mu_sigma.LOS = 6.8;
paratable(4).mu_sigma.MD  = 9.4;
paratable(4).mu_sigma.SD  = 6.3;

paratable(4).mu_c.LOS = 7.2;
paratable(4).mu_c.MD  = 5.4;
paratable(4).mu_c.SD  = 4.9;

paratable(4).d_c_min.LOS = 4.4;
paratable(4).d_c_min.MD  = 1.1;
paratable(4).d_c_min.SD  = 1.0;

paratable(4).y1.SD = -10;
paratable(4).y1.DI = -13.5;

paratable(4).y2.SD = 10;
paratable(4).y2.DI = 13.5;

paratable(4).W_DI    = 5;
paratable(4).W_road  = 18;
paratable(4).N_lanes = 4;