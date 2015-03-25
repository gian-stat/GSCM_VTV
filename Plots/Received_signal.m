function [ graph ] = Received_signal( s, y, ts )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
graph=s;
t=0:size(y,1)-1;
N_tx=size(y,2);
N_rx=size(y,3);
nfig=2+3*N_tx+N_rx;
for w=1:N_tx 
    for n=1:N_rx
        figure(nfig);
        subplot (3,1,1), stem(t*ts*1e6,(abs(s)))
        title(['Received signal - Channel ','TX',num2str(w),'RX',num2str(n)]);
        xlabel('Time [micros]');
        ylabel('Transmitted signal');
        hold on
        grid on
        
        subplot (3,1,2), stem(t*ts*1e6,(abs(y(:,w,n))))
        xlabel('Time [micros]');
        ylabel('Received signal amplitude');
        hold on
        grid on
        subplot (3,1,3), stem(t*ts*1e6,((angle(y(:,w,n)))))
        
        xlabel('Time [micros]');
        ylabel('Received signal phase');
        hold on
        grid on
        fixfig;
        nfig=nfig+1;
    end
end

