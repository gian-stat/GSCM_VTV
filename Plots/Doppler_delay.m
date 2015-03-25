function [ doppler ] = Doppler_delay( h_TOT,d_ni,ni_max, d_tau)
%DOPPLER DELAY plot Doppler shift vs delay
%   Detailed explanation goes here
doppler=h_TOT;
h_DD=permute(sum(h_TOT,1),[3 2 4 5 1]);
dop=-ni_max:d_ni:ni_max;
tau=(0:size(h_TOT,2))*d_tau*1e6;
for w=1:size(h_TOT,4)
    figure(1+2*size(h_TOT,4)+size(h_TOT,5)+w)
    for v=1:size(h_TOT,5)  
        subplot(size(h_TOT,5),1,v),imagesc (tau,dop,(20*log10(abs(h_DD(:,:,w,v)))))
        
        hold on
        grid on
        ylabel('Doppler shift [Hz]');
        xlabel('Delay [µs]');
        title(['Delay Doppler Spectrum - ','TX',num2str(w),' - RX',num2str(v)]);
        colorbar;
        colormap(hot);
        fixfig;
    end
end

end

