
%% %%%%%%%%%%%%%%%%%%%%%%% MASTER THESIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A MATLAB implementation of a Geometry Based Stochastic Channel Model %%%
%% %%%%%%%%%%%%% for Vehicle to Vehicle communication %%%%%%%%%%%%%%%%%%%%%

% Student: Gianluca Statuti
% Examiner: Prof. Erik Strom
% Supervisor: Keerthi Nagalapur
%% References
% [1] Gianluca Statuti, "A MATLAB implementation of a geometry-based
% stochastic channel model", Master Thesis, 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE: it is necessary to include all folders in MATLAB path

%% %%%%%%%%%%%%%%%%%%% START of IMPLEMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all                             % close all graphic windows
clear                                 % reset all variables
clc                                   % clear screen
format long                           % maximum precision in representing 
                                      % numbers
tic
%% %%%%%%%%%%%%%%%%%%%%%%% INPUT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] - (Appendix A)

% General parameters
f_0=5.9e9;               % central frequence (Hz) - [1]
c_0=3e8;                 % light speed (m/s)
kw= 2*pi*f_0/c_0;        % wave number

%%%%%%%%%%%%%%%%%%%%%%%%%% ANTENNA SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING ANTENNA PATTERN
% Please select a .mat file containing the pattern vector
pattern_p=importdata('Isotropic_2D_pattern.mat');
% we could use different antenna patterns for Tx and Rx importing different
% patterns and changing the assignment below

% Resampling antenna pattern
r_scale=1; % scaling factor: i.e. new_pattern_vector_length=ceil(old_pattern_vector_length/r_scale)
pattern_tx_a=downsample(pattern_p,r_scale); % Transmitter antenna pattern
pattern_rx_a=downsample(pattern_p,r_scale); % Receiver antenna pattern

% MIMO PARAMETERS
N_a_tx=1;          % number of antennas of TX
N_a_rx=1;          % number of antennas of RX
delta_y_tx=c_0/f_0*0.5;    % distance between antennas elements at TX
delta_y_rx=c_0/f_0*0.5;    % distance between antennas elements at RX

pattern_tx{N_a_tx}=[];
for w=1:N_a_tx
    pattern_tx{w}=pattern_tx_a;
end
pattern_rx{N_a_rx}=[];
for v=1:N_a_rx
    pattern_rx{v}=pattern_rx_a;
end

% Antenna gain
g_T_dB=[3 3];        % TX antenna gain [dB]
g_R_dB=[3 3];        % RX antenna gain [dB]

g_T=10.^(g_T_dB/10); % TX antenna gain
g_R=10.^(g_R_dB/10); % RX antenna gain

% Uncomment the following line to plot TX/RX antenna pattern
% pattern_plot(pattern_tx_a,pattern_rx_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% GEOMETRY SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lenght of road strip (m)
l_road = 200;    % input

% Select the kind of environment: %%%%%%%%%%%%%%

environment = 0; % Rural=0 Highway=1

% Select the kind of TX/RX positions assignment: %%%%%%%%%%%%%%

position = 0;    % Manual=0 Statistical=1

% TX/RX initial coordinates
% NOTE: it is valid only if position=0
x_rm = 170;  % [m] [0, l_road]
x_tm = 95;   % [m] [0, l_road]
y_rm = 2;    % [m] e.g. n*W_road/(N_lanes*2) with n integer in [-N_lanes/2, N_lanes/2]
y_tm = -2;   % [m] e.g. n*W_road/(N_lanes*2) with n integer in [-N_lanes/2, N_lanes/2]

% TX/RX speeds
% NOTE: use absolute values
v_t = 30;    % [m/s]
v_r = 40;    % [m/s]

% Select the direction of the travel: %%%%%%%%%%%%%%
% NOTE: it is valid only if we decide to assign TX and RX initial position
% in a stochastic way

direction = 1;     % SM=0 OP=1

% Select kind of simulation

simulation = 1;   % 0 = no input signal 1 = input signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% TIME, DELAY and DOPPLER SHIFT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
% Doppler shift
N_d = 100;                     % (input) number of Doppler shift quantization level

v_max=max(abs([v_r,v_t]));    % maximum speed [m/s]
ni_max=((2*v_max)/c_0)*f_0;   % maximum doppler shift quantized [Hz]
d_ni=2*ni_max/N_d;            % width of Doppler shift bin [s] [Hz]

Tc=1/(2*ni_max);              % coherence time [s] - [1] - eq. (3.11)

if simulation==0              % [1] - eq. (3.12)
    Ts=Tc/2;                  % h sample time [s] (=Tc)
    f_l=10*Ts;               % frame length [s]
    d_tau=0.01e-7;            % (input) width of delay bin [s]
    N_samp=floor(f_l/Ts);     % number of samples
else
    
    %%%%%%%%%%%%%%%%%%%%%%%% INPUT SIGNAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
    ts=(0.1e-6)/2;           % signal sample time [s] (<=d_tau) - [1] - eq. (4.2)
    s_l=10e-4;               % signal length [s]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delay resolution
    % NOTE: typical max delay is around 0.5-1 µs depending on the length of
    % road strip we consider - [1] - eq. (4.5)
    d_tau=ts;                     % (input) width of delay bin [s] (>=ts)
    
    % Time settings
    % Channel's characteristics change significantly only every ~Tc (coherence
    % time). For that reason we sample the channel at the first sampling rate
    % that is a multiple integer of ts (signal sampling time) but less than
    % equal of an arbitrary fraction of Tc.
    
    rf=100;                      % (input) reduction factor for Tc
    
    % Channel's sampling rate selection - [1] - eq. (4.3)
    uf=1;                       % upsampling factor initialization
    Ts=ts;                      % Ts initialization
    while Ts<Tc/rf              % make Ts<Tc/rf but integer multiple of ts
        Ts=uf*ts;               % h sample time [s]
        uf=uf+1;
    end
    %%%%%%%
    f_l=s_l;                    % frame length [s]
    N_samp=floor(f_l/Ts);       % number of samples
    N_samp_sig=(uf-1)*N_samp;   % number of samples
end

% Other parameters
f_s= 1/Ts;                  % channel's sampling frequence (Hz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF INPUT SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%% PARAMETERS ASSIGNMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Channel's parameters are imported from the table - [1] - Table 2.1
Table1;                   % build parameters table
it = environment+1;       % paratable index

% Road characteristics

N_lanes = paratable(it).N_lanes;    % number of lanes
W_road  = paratable(it).W_road;     % road width
W_DI    = paratable(it).W_DI;        % diffuse scatterer width

% Number of scatterers
P = ceil(paratable(it).chi.MD*l_road); % number of MD scatterers
Q = ceil(paratable(it).chi.SD*l_road); % number of SD scatterers
R = ceil(paratable(it).chi.DI*l_road); % number of DI scatterers

% Scatterer means positions

y_1SD = paratable(it).y1.SD; % means y-coordinates of SD scatterers
y_2SD = paratable(it).y2.SD; % means y-coordinates of SD scatterers

y_1DI = paratable(it).y1.DI; % means y-coordinates of DI scatterers
y_2DI = paratable(it).y2.DI; % means y-coordinates of DI scatterers

% Other parameters

n_MD  = paratable(it).n.MD;          % pathloss exponent vector MD
n_SD  = paratable(it).n.SD;          % pathloss exponent vector SD
n_LOS = paratable(it).n.LOS;         % pathloss exponent vector LOS
n_DI  = paratable(it).n.DI;          % pathloss exponent vector DI

mu_sigma = table2array(struct2table(paratable(it).mu_sigma)); % lowest value variance vector LOS MD SD
mu_c     = table2array(struct2table(paratable(it).mu_c));     % lowest value distance vector LOS MD SD

d_cmin   = table2array(struct2table(paratable(it).d_c_min));  % minimum 0.5-coherence distance vector LOS MD SD

G_0_LOS_dB = paratable(it).G0.LOS;       % reference power vector LOS [dB]
G_0_SD_dB  = paratable(it).G0.SD;        % reference power vector SD [dB]
G_0_MD_dB  = paratable(it).G0.MD;        % reference power vector MD [dB]
G_0_DI_dB  = paratable(it).G0.DI;        % reference power vector DI [dB]

G_0_LOS=10^(G_0_LOS_dB/10);  % reference power vector LOS
G_0_SD=10.^(G_0_SD_dB/10);   % reference power vector SD
G_0_MD=10.^(G_0_MD_dB/10);   % reference power vector MD
G_0_DI=10^(G_0_DI_dB/10);    % reference power vector DI

% Maximum delay and number of delay's quantization level %%%%%%%%%%%%%%%%%%
tau_max=(2*(sqrt((l_road).^2+(y_2DI+W_DI/2).^2)))/c_0; % maximum quantized delay [s]
N = floor(tau_max/d_tau);   % number of delay's quantization level

%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%% % INITIAL POSITION/SPEEDS OF TX, RX, SD, MD, DI SCATTERERS %%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% TX/RX INITIAL COORDINATES %%%%%%%%%%%%%%%%%%%
if position==0 % manual
    % MANUAL ASSIGNMENT
    x_r0=x_rm;
    x_t0=x_tm;
    y_r0=y_rm;
    y_t0=y_tm;
else % statistical
    % STATISTICAL ASSIGNMENT
    % Tx and Rx starting coordinates are assigned in a statistical way
    x_r0=rand(1,1)*l_road;
    x_t0=rand(1,1)*l_road;
    
    % 2 random lane's index
    y_plus=2*floor(rand(1,1)*N_lanes/2)+1; % NOTE: the lane is stochastically selected
    y_minus=-(2*floor(rand(1,1)*N_lanes/2)+1);
    
    w=rand(1,1);    % random value used to define Tx and Rx way
    if direction==0
        if w>0.5
            y_r0=abs(y_plus*W_road/(2*N_lanes)); % y-coordinates of Tx and Rx has the same positive sign
            y_t0=abs(y_minus*W_road/(2*N_lanes));
        else
            y_r0=-abs(y_plus*W_road/(2*N_lanes)); % y-coordinates of Tx and Rx has the same negative sign
            y_t0=-abs(y_minus*W_road/(2*N_lanes));
        end
    else
        if w>0.5
            y_r0=y_plus*W_road/(2*N_lanes); % y-coordinates of Tx and Rx has different signs
            y_t0=y_minus*W_road/(2*N_lanes); % Rx positive - Tx Negative
        else
            y_r0=-y_plus*W_road/(2*N_lanes); % y-coordinates of Tx and Rx has different signs
            y_t0=-y_minus*W_road/(2*N_lanes); % Tx positive - Rx Negative
        end
    end
end

% COORDINATES OF MIMO ANTENNAS ARRAY
% We are considering only linear antennas arrays orthogonal relative to the
% direction of the travel

% TX ANTENNAS ARRAY - [1] - eq. (3.18)
y_at=y_t0+(-(N_a_tx-1)*(delta_y_tx/2):delta_y_tx:(N_a_tx-1)*(delta_y_tx/2));

% RX ANTENNAS ARRAY - [1] - eq. (3.19)
y_ar=y_r0+(-(N_a_rx-1)*(delta_y_rx/2):delta_y_rx:(N_a_rx-1)*(delta_y_rx/2));

% Initial LOS distance computed from the array's center
d_LOS0=sqrt((abs(x_r0-x_t0)).^2 + (abs(y_r0-y_t0).^2));

%%%%%%%%%%%%%%%%%%% SCATTERERS DISTRIBUTION %%%%%%%%%%%%%%%%%%%

% MD scatterers %%%%%%%%%%%%%%%% [1] - Section 2.2.1

% Initials coordinates generation for P MD scatterers
x_MD_0=rand(1,P)*l_road;
y_MD_0=zeros(1,ceil(P));

y_p=2*floor(rand(1,P)*N_lanes/2)+1; % NOTE: the lane is stochastically selected

z=rand(1,P);
for j=1:(ceil(P))
    if z(j)>0.5
        y_MD_0(j)=abs(y_p(j)*W_road/(2*N_lanes)); % y-coordinates of MD
    else
        y_MD_0(j)=-abs(y_p(j)*W_road/(2*N_lanes));
    end
end

% SD scatterers %%%%%%%%%%%%%%%% [1] - Section 2.2.1

% initialization of coordinates and standard deviation

sigma_ySD=sqrt((y_2DI-(W_DI/2)-W_road/2)/2); % [1] - eq. (3.1)

% Coordinates generation for Q SD scatterers
x_SD=rand(1,Q)*l_road;

y_SD(1:ceil(Q/2))=y_1SD+randn(1,ceil(Q/2))*sigma_ySD;
y_SD(ceil(Q/2)+1:Q)=y_2SD+randn(1,floor(Q/2))*sigma_ySD;

% DI scatterers %%%%%%%%%%%%%% [1] - Section 2.2.1

x_DI=rand(1,R)*l_road;

y_DI(1:ceil(R/2))=(y_1DI-W_DI/2)+rand(1,ceil(R/2))*W_DI;
y_DI(ceil(R/2)+1:R)=(y_2DI-W_DI/2)+rand(1,floor(R/2))*W_DI;

%%% SPEED VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tx and Rx absolute speeds
% NOTE: We are assuming a rigth-hand traffic
signRX=-sign(y_r0);         % direction of speed vector
v_xr=signRX.*v_r;           % absolute speeds

signTX=-sign(y_t0);         % direction of speed vector
v_xt=signTX.*v_t;           % absolute speeds

% MD scatterers speeds [1] - eq. (3.2 - 3.3 - 3.4)
% NOTE: parametrize mean value and variance
if environment==0 % rural
    v_MD=abs(20+randn (1,P)*5);  % stochastic generated speeds
else % highway
    v_MD=abs(30+randn (1,P)*5);  % stochastic generated speeds
end
signMD=-sign(y_MD_0);        % direction of speed vector
v_xMD=signMD.*v_MD;          % absolute speeds

%% %%%%%%%%%%%%%% INITIALS VIDEO PRINTS %%%%%%%%%%%%%%%%
if environment ==0
    fprintf('%s\n','Kind of environment: Rural')
else    fprintf('%s\n','Kind of environment: Highway')
end

fprintf('%s%d\n','Number of lanes: ',N_lanes)
fprintf('%s%d%s\n','Road length: ',l_road,' m')
fprintf('%s%d%s\n','Road width: ',W_road,' m')
fprintf('%s%d\n','Number of SD scatterers: ',Q)
fprintf('%s%d\n','Number of MD scatterers: ',P)
fprintf('%s%d\n','Number of DI scatterers: ',R)

fprintf('%s%d%s\n','Frame length: ',f_l*10^3,' ms')
fprintf('%s%d%s\n','Sample time: ',Ts*10^6,' µs')

fprintf('%s%d\n','Number of TX antennas elements: ', N_a_tx)
fprintf('%s%d\n','Number of RX antennas elements: ', N_a_rx)

fprintf('%s%d%s\n','Average initial LOS distance: ',d_LOS0,' m')

fprintf('%s%d%s\n','TX speed: ',v_xt*3.6,' km/h')
fprintf('%s%d%s\n','RX speed: ',v_xr*3.6,' km/h')


%% Variables initialization

x_r=zeros(1,N_samp);
x_t=zeros(1,N_samp);
y_r=zeros(1,N_samp);
y_t=zeros(1,N_samp);

x_MD=zeros(P,N_samp);
y_MD=zeros(P,N_samp);

AOD_LOS=zeros(1,N_samp);
AOA_LOS=zeros(1,N_samp);
AOD_LOS_p=zeros(1,N_samp);
AOA_LOS_p=zeros(1,N_samp);
AOD_LOS_q=zeros(1,N_samp);
AOA_LOS_q=zeros(1,N_samp);

a_LOS=zeros(1,N_samp);
a_SD=zeros(Q,N_samp);
a_MD=zeros(P,N_samp);
a_DI=zeros(R,N_samp);

d_LOS=zeros(1,N_samp);
tau_LOS=zeros(1,N_samp);
tau_SD=zeros(Q,N_samp);
tau_MD=zeros(P,N_samp);
tau_DI=zeros(R,N_samp);

tau_LOS_q=zeros(1,N_samp);

ni_LOS=zeros(1,N_samp);
ni_SD=zeros(Q,N_samp);
ni_MD=zeros(P,N_samp);
ni_DI=zeros(R,N_samp);

ni_LOS_q=zeros(1,N_samp);

h_LOS=zeros(N_samp,N+1,N_d+1,N_a_tx,N_a_rx);
h_SD=zeros(N_samp,N+1,N_d+1,N_a_tx,N_a_rx);
h_MD=zeros(N_samp,N+1,N_d+1,N_a_tx,N_a_rx);
h_DI=zeros(N_samp,N+1,N_d+1,N_a_tx,N_a_rx);

d_ref=1;

%% COMPLEX PARAMETERS COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VARIANCE of AMPLITUDE GAIN %%% [1] - Section 2.2.3
sigma_S_LOS=exprnd(mu_sigma(1,1));          % variance for LOS

sigma_S_MD=exprnd(mu_sigma(1,2),[1 P]);     % variance for MD

sigma_S_SD=exprnd(mu_sigma(1,3),[1 Q]);     % variance for SD

sigma_cr=1;                                 % variance for DI 

% COHERENCE DISTANCE %%% [1] - Section 2.2.3
d_c_LOS=d_cmin(1,1)+exprnd(mu_c(1,1));      % coherence distance for LOS

d_c_MD=d_cmin(1,2)+exprnd(mu_c(1,2),[1 P]); % coherence distance for MD

d_c_SD=d_cmin(1,3)+exprnd(mu_c(1,3),[1 Q]); % coherence distance for SD

% AMPLITUDE COMPUTING PARAMETERS %%%

% PATH GAIN %%% [1] - eq (3.7) - Appendix A
g_S_LOS=gain_computation (sigma_S_LOS,d_c_LOS); % LOS path gain

g_S_MD= zeros (1,P);
for j=1:P
    g_S_MD(j)=gain_computation (sigma_S_MD(j),d_c_MD(j)); % MD path gain
end

g_S_SD= zeros (1,Q);
for j=1:Q
    g_S_SD(j)=gain_computation (sigma_S_SD(j),d_c_SD(j)); % SD path gain
end

% DI scatterers gain - [1] - Section 2.2.2
c_r=sqrt(1/2)*randn(1,R)+1i*sqrt(1/2)*randn(1,R);

% PHASEs %% [1] - Section 2.2.2
phi_LOS= rand (1,1)*(2*pi); % LOS

phi_SD=rand (1,Q)*(2*pi);   % SD SCATTERERS

phi_MD= rand (1,P)*(2*pi);  % MD SCATTERERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% START OF MIMO CYCLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for v=1:N_a_rx % RX ANTENNAS ELEMENT CYCLE
    y_r0=y_ar(v); % set RX antenna initial position
    for w=1:N_a_tx % TX ANTENNAS ELEMENT CYCLE
        y_t0=y_at(w); % set TX antenna initial position
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%%%%%%%%%%%%%%%%%% START OF SAMPLING CYCLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i=1:N_samp; % SAMPLING CYCLE
            %%
            %%%% TIME VARYING POSITION OF TX, RX, SD, MD, DI SCATTERERS %%%
            
            %%%% TX and RX %%%%%
            
            % Time-varying positions of Tx, Rx and scatterers
            x_r(i)=x_r0 + (v_xr) *i*Ts;
            x_t(i)=x_t0 + (v_xt) *i*Ts;
            y_r(i)=y_r0;
            y_t(i)=y_t0;
            
            % Time-varying positions of MD scatterers
            x_MD(:,i)=x_MD_0 + v_xMD*i*Ts;
            y_MD(:,i)=y_MD_0;
          
            %%
            %%%%%% DISTANCE, AOD and AOA COMPUTATION %%%%%%%%%%%%%%%%%%%%%
            
            %%% LOS %%%
            % DISTANCE - [1] - eq. (3.13)
            d_LOS(i)=sqrt((abs(x_r(i)-x_t(i))).^2 + (abs(y_r(i)-y_t(i)).^2)); % distance of LOS path
           
            % [1] - Appendix A
            % AOD (rad) - Computed starting from TX speed vector direction
            [AOD_LOS(i),AOD_LOS_p(i)]= AOD_computation (x_t(i), y_t(i), x_r(i), y_r(i),v_xt);
            AOD_LOS_q(i)=Angle_index(AOD_LOS(i),length(pattern_tx{w}));
            % AOA (rad) - Computed starting from RX speed vector direction
            [AOA_LOS(i),AOA_LOS_p(i)]= AOA_computation (x_t(i), y_t(i), x_r(i), y_r(i),v_xr);
            AOA_LOS_q(i)=Angle_index(AOA_LOS(i),length(pattern_rx{v}));
            
            %%% MD scatterers %%%
            AOD_MD=zeros(1,P);      % MD AOD initialization
            AOA_MD=zeros(1,P);      % MD AOA initialization
            AOD_MD_p=zeros(1,P);    % MD AOD initialization
            AOA_MD_p=zeros(1,P);    % MD AOD initialization
            
            % MD distance and angles computations - [1] - eq. (3.13)
            d_TX_MD=sqrt((abs(x_t(i)-x_MD(:,i))).^2 + (abs(y_t(i)-y_MD(:,i)).^2)); % TX -> MD path
            d_MD_RX=sqrt((abs(x_r(i)-x_MD(:,i))).^2 + (abs(y_r(i)-y_MD(:,i)).^2)); % MD -> RX path
            
            % [1] - Appendix A
            for j=1:ceil(P)
                [AOD_MD(j),AOD_MD_p(j)]= AOD_computation (x_t(i), y_t(i), x_MD(j,i), y_MD(j,i),v_xt);% (rad) calculated starting from speed vector direction
                [AOA_MD(j),AOA_MD_p(j)]= AOA_computation (x_MD(j,i), y_MD(j,i), x_r(i), y_r(i),v_xr);% (rad) calculated starting from speed vector direction
            end
            AOD_MD_q=Angle_index(AOD_MD,length(pattern_tx{w})); % quantized AOD angle
            AOA_MD_q=Angle_index(AOA_MD,length(pattern_rx{v})); % quantized AOA angle
            
            %%% SD scatterers %%%
            AOD_SD=zeros(1,Q);      % SD AOD initialization
            AOA_SD=zeros(1,Q);      % SD AOD initialization
            AOD_SD_p=zeros(1,Q);    % SD AOD initialization
            AOA_SD_p=zeros(1,Q);    % SD AOD initialization
            
            % SD distance and angles computations - [1] - eq. (3.13)
            d_TX_SD=sqrt((abs(x_t(i)-x_SD)).^2 + (abs(y_t(i)-y_SD).^2)); % TX -> SD path
            d_SD_RX=sqrt((abs(x_r(i)-x_SD)).^2 + (abs(y_r(i)-y_SD).^2)); % SD -> RX path
            
            % [1] - Appendix A
            for j=1:ceil(Q)
                [AOD_SD(j),AOD_SD_p(j)]= AOD_computation (x_t(i), y_t(i), x_SD(j), y_SD(j),v_xt);% (rad) calculated starting from speed vector direction
                [AOA_SD(j),AOA_SD_p(j)]= AOA_computation (x_SD(j), y_SD(j), x_r(i), y_r(i),v_xr);% (rad) calculated starting from speed vector direction
            end
            
            AOD_SD_q=Angle_index(AOD_SD,length(pattern_tx{w})); % quantized AOD angle
            AOA_SD_q=Angle_index(AOA_SD,length(pattern_rx{v})); % quantized AOA angle
            
            %%% DI scatterers %%%
            AOD_DI=zeros(1,R);      % DI AOD initialization
            AOA_DI=zeros(1,R);      % DI AOA initialization
            AOD_DI_p=zeros(1,R);    % DI AOD initialization
            AOA_DI_p=zeros(1,R);    % DI AOD initialization
            
            % DI scatterer on the lower part of the road - [1] - eq. (3.13)
            d_TX_DI=sqrt((abs(x_t(i)-x_DI)).^2 + (abs(y_t(i)-y_DI).^2)); % TX -> DI path
            d_DI_RX=sqrt((abs(x_r(i)-x_DI)).^2 + (abs(y_r(i)-y_DI).^2)); % DI -> RX path
            
            % [1] - Appendix A
            for j=1:ceil(R)
                [AOD_DI(j),AOD_DI_p(j)]= AOD_computation (x_t(i), y_t(i), x_DI(j), y_DI(j),v_xt);% (rad) calculated starting from speed vector direction
                [AOA_DI(j),AOA_DI_p(j)]= AOA_computation (x_DI(j), y_DI(j), x_r(i), y_r(i),v_xr);% (rad) calculated starting from speed vector direction
            end
            AOD_DI_q=Angle_index(AOD_DI,length(pattern_tx{w})); % quantized AOD angle
            AOA_DI_q=Angle_index(AOA_DI,length(pattern_rx{v})); % quantized AOA angle
            
            %%
            
            %%%%%%%%%% AMPLITUDE COMPUTATION %%% [1] - eq. (3.5) %%%%%%%
            a_LOS(i)=g_S_LOS*exp(1i*phi_LOS)*G_0_LOS^(1/2)*(d_ref/d_LOS(i))^(n_LOS/2); % LOS amplitude vector
            
            a_MD(:,i)=g_S_MD.*exp(1i*phi_MD).*G_0_MD.^(1/2).*(d_ref/(d_TX_MD+d_MD_RX)).^(n_MD/2); % MD amplitude vector
            
            a_SD(:,i)=g_S_SD.*exp(1i*phi_SD).*G_0_SD.^(1/2).*(d_ref./(d_TX_SD+d_SD_RX)).^(n_SD/2); % SD amplitude vector
            
            % [1] - eq. (3.8)
            a_DI(:,i)=G_0_DI^(1/2)*c_r.*(d_ref./(d_TX_DI.*d_DI_RX)).^(n_DI/2); % DI amplitude vector
           
            %%
            %%%%%%%%%%%%%% DELAY COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Propagation delays - [1] - eq. (3.14) (Appendix A)
            tau_LOS(i)=(d_LOS(i)/c_0);                  % LOS delay
            tau_LOS_q(i)=Delay_index(tau_LOS(i),N,tau_max);   % quantized LOS delay
            
            tau_MD(:,i)=((d_TX_MD+d_MD_RX)./c_0);       % MD delay
            tau_MD_q=Delay_index(tau_MD(:,i),N,tau_max);      % quantized MD delay
            
            tau_SD(:,i)=((d_TX_SD+d_SD_RX)./c_0);       % SD delay
            tau_SD_q=Delay_index(tau_SD(:,i),N,tau_max);      % quantized SD delay
            
            tau_DI(:,i)=((d_TX_DI+d_DI_RX)./c_0);       % DI delay
            tau_DI_q=Delay_index(tau_DI(:,i),N,tau_max);      % quantized DI delay
            
            %%
            %%%%%%%%%%%%%% DOPPLER SHIFTS COMPUTATION %%%%%%%%%%%%%%%%%%%%%
            % [1] - eq. (3.15) (Appendix A)
            
            ni_LOS(i)=f_0/c_0*((v_xt-v_xr)*cos(AOD_LOS_p(i)));                            % LOS Doppler shift
            ni_LOS_q(i)=Doppler_index(ni_LOS(i),N_d,ni_max);                              % quantized LOS Doppler shift
            
            ni_MD(:,i)=f_0/c_0*((v_xt-v_xMD).*cos(AOD_MD_p)+(v_xr-v_xMD).*cos(AOA_MD_p)); % MD Doppler shift
            ni_MD_q=Doppler_index(ni_MD(:,i),N_d,ni_max);                                 % quantized MD Doppler shift
            
            ni_SD(:,i)=f_0/c_0*(v_xt*cos(AOD_SD_p)+v_xr*cos(AOA_SD_p));                   % SD Doppler shift
            ni_SD_q=Doppler_index(ni_SD(:,i),N_d,ni_max);                                 % quantized SD Doppler shift
            
            ni_DI(:,i)=f_0/c_0*(v_xt*cos(AOD_DI_p)+v_xr*cos(AOA_DI_p));                   % DI Doppler shift
            ni_DI_q=Doppler_index(ni_DI(:,i),N_d,ni_max);                                 % quantized DI Doppler shift
            
            %%
            %%% TIME-DELAY_DOPPLER TRANSFER FUNCTION COMPUTATION %%%%%%%%%%
            % [1] - eq. (3.16)
            
            h_LOS(i,tau_LOS_q(i),ni_LOS_q(i),w,v)=a_LOS(i)*exp(1i*kw*d_LOS(i))*g_R(v)*pattern_rx{v}(AOA_LOS_q(i))*g_T(w)*pattern_tx{w}(AOD_LOS_q(i));
            for j=1:Q
                h_SD(i,tau_SD_q(j),ni_SD_q(j),w,v)=h_SD(i,tau_SD_q(j),ni_SD_q(j),w,v)+a_SD(j)*exp(1i*kw*(d_TX_SD(j)+d_SD_RX(j)))*g_R(v)*pattern_rx{v}(AOA_SD_q(j))*g_T(w)*pattern_tx{w}(AOD_SD_q(j));
            end
            
            for j=1:P
                h_MD(i,tau_MD_q(j),ni_MD_q(j),w,v)=h_MD(i,tau_MD_q(j),ni_MD_q(j),w,v)+a_MD(j)*exp(1i*kw*(d_TX_MD(j)+d_MD_RX(j)))*g_R(v)*pattern_rx{v}(AOA_MD_q(j))*g_T(w)*pattern_tx{w}(AOD_MD_q(j));
            end
            for j=1:R
                h_DI(i,tau_DI_q(j),ni_DI_q(j),w,v)=h_DI(i,tau_DI_q(j),ni_DI_q(j),w,v)+a_DI(j)*exp(1i*kw*(d_TX_DI(j)+d_DI_RX(j)))*g_R(v)*pattern_rx{v}(AOA_DI_q(j))*g_T(w)*pattern_tx{w}(AOD_DI_q(j));
            end
            
        end % end of time cycle
        
        %% %%%%%%%%%%% GEOMETRY PLOT %%%%%%%%%%%%%%% [1] - Appendix A
        % Uncomment to plot the geometry of the propagation environment
        geometry(x_t,x_r,y_t,y_r,x_SD,y_SD,x_MD,y_MD,x_DI,y_DI,Q,R, l_road,W_road, environment);
    
    end % end of TX element cycle
end % end of RX element cycle

% TOTAL TIME-DELAY_DOPPLER TRANSFER FUNCTION COMPUTATION %%%%%%%%%%
% [1] - eq. (3.17)
h_TOT=h_LOS+h_SD+h_MD+h_DI;

h_TOT_p=permute(sum(h_TOT,3), [1 2 4 5 3]); % sum up all Doppler shift components
% END OF h(t,tau) COMPUTATION

%% %%%%%%%%%%%%%%%%%% TRANSCEIVER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if simulation==1
    %%%%%%%%%%%%%%%%%%%%% TEST SIGNAL GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Impulse train - [1] - eq. (4.6)
    T=N+1; % period of impulse train (must be >tau_max: i.e. T=2*tau_max)
    
    s_in=zeros(1,N_samp_sig);% generate an impulse train of length N_samp_sig
    s_in(1:T:N_samp_sig)=1;  % and period T
    
    % Delayed signal
    % Build a matrix of the input signal translated and padded with 0
    s=zeros(N_samp_sig,N+1);
    for j=0:N
        s(:,j+1)=padarray(s_in(1:end-j),[0 j],0,'pre');
    end
    %%%%%%%%%%%%%%%%%% UPSAMPLE and OUTPUT COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%
    h=zeros(N_samp_sig,N+1,N_a_tx,N_a_rx);      % initialization
    o=ones (N_samp_sig,N+1);             % initialization
    y=zeros(N_samp_sig,N_a_tx,N_a_rx);          % initialization
    
    for w=1:N_a_tx % for each TX antenna element
        for v=1:N_a_rx % for each RX antenna element
            for k=1:N+1 % for each delay bin
                h(:,k,w,v)=interp(h_TOT_p(:,k,w,v),uf-1); % upsampling h(t,tau)
            end
            
            % Filtering %%%%%%%%%%% - [1] - eq. (4.6)
            for k=1:N+1
                y(:,w,v)= y(:,w,v)+s(:,k).*h(:,k,w,v);
            end
        end
    end

    %% %%%%%%%%%%%% OUTPUT PLOT %%%%%%%%%%%%%%%%% [1] - Appendix A
    % Uncomment to plot transmitted signal and abs/angle of received signal from each
    % antenna element 
    Received_signal(s_in,y,ts);
    %%%%%%%%%%%%%
end
%% %%%%%%%%%%% h(t,tau) PLOTS %%%%%%%%%%%%%%% [1] - Appendix A
% % Uncomment to plot time-delay impulse response for each RX antenna element
% Time_delay(h_TOT_p,Ts,d_tau);
% 
% % Uncomment to plot Doppler-time impulse response for each RX antenna element
% Doppler_time(h_TOT,Ts,d_ni,ni_max);
% 
% % Uncomment to plot Doppler-delay impulse response for each RX antenna element
% Doppler_delay(h_TOT,d_ni,ni_max,d_tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display h(t,tau) evolution in time
global h_T;
global t_sample;
global d_t;
h_T=h_TOT_p(:,:,1,1); % first TX and first RX element
t_sample = Ts;
d_t = d_tau;
% Uncomment to open GUI
Impulse_res;
%%%%%%%%%%%%%
toc