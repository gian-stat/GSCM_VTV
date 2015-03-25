function [ g_s ] = gain_computation( sigma_s, d_c )
% This function provide a random value of g_s given by a distribution
% obtained correlating a normal distributed random sequence with a Gaussian
% autocorrelation function by means of convolution
N=250;                  % number of samples of distributions
v=randn(1,N);           % generation of normal distributed random vector          

delta_d=0:0.1:25; % generation of the set of values in which we calculate the autocorrelation function
r_d=sigma_s*exp(-((log(2))/(d_c^2)).*(delta_d).^2); % autocorrelation function

% Convolution
r_d_n=abs(r_d)/abs(r_d(1)); % normalize r_d
G_s_dB=conv(r_d_n,v);        % convolution between the two distribution

% CONVERT IN LINEAR
G_s=10.^(G_s_dB/10);

g_s=G_s(randi(N,1));    % pick one of the previously generated samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

