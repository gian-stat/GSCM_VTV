function [ delay ] = Time_delay( h_TOT_p,Ts,d_tau)
% TIME DELAY plots the time-delay impulse response
% For each RX antenna element abs and angle of h(t,tau) are plotted

delay=h_TOT_p;
tau=(0:size(h_TOT_p,2))*d_tau*1e6;
t=(0:size(h_TOT_p,1))*Ts;
nfig=2;
for w=1:size(h_TOT_p,3)
    for v=1:size(h_TOT_p,4)
        figure(nfig);
        subplot(2,1,1), imagesc (tau,t,(20*log10(abs(h_TOT_p(:,:,w,v)))))
        hold on
        grid on
        colorbar;
        colormap(hot);
        xlabel('Propagation delay [µs]');
        ylabel('Time [s]');
        title({['Channel ', 'TX',num2str(w),' RX', num2str(v)];'|h(t,\tau)|^2 [dB]'});
        
        subplot(2,1,2), imagesc (tau,t,(angle(h_TOT_p(:,:,w,v))))
        hold on
        grid on
        colorbar;
        colormap(hot);
        xlabel('Propagation delay [µs]');
        ylabel('Time [s]');
        title('Phase of h(t,\tau)');
        fixfig;
        nfig=nfig+1;
    end
end
end

