function [ plot ] = pattern_plot( pattern_tx, pattern_rx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
plot=pattern_tx;
theta_tx=0:2*pi/length(pattern_tx):2*pi;
theta_rx=0:2*pi/length(pattern_rx):2*pi;

figure(20);
subplot(1,2,1), polar (theta_tx(1:end-1),pattern_tx)
hold on
grid on
title('TX antenna pattern')
subplot(1,2,2), polar (theta_rx(1:end-1),pattern_rx)
hold on
grid on
title('RX antenna pattern')
end

