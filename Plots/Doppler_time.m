function [ doppler ] = Doppler_time( h_TOT,Ts,d_ni,ni_max)
%TIME DELAY plot time vs Doppler shift
doppler=h_TOT;
h_DT=permute(sum(h_TOT,2),[1 3 4 5 2]);

t=(0:size(h_TOT,1))*Ts;
dop=-ni_max:d_ni:ni_max;
for w=1:size(h_TOT,4)
    figure(1+size(h_TOT,4)+size(h_TOT,5)+w)
    for v=1:size(h_TOT,5)
        
        subplot(size(h_TOT,5),1,v),imagesc (dop,t,(20*log10(abs(h_DT(:,:,w,v)))))
        
        hold on
        grid on
        xlabel('Doppler shift [Hz]');
        ylabel('Time [s]');
        title(['Doppler Power Spectral Density - ','TX',num2str(w),' - RX',num2str(v)]);
        colorbar;
        colormap(hot);
    end
    fixfig;
end
end

