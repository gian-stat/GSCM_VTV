function [ ind ] = Delay_index( tau, N, tau_max )
% index compute the level's index associated to a tau
d_tau=tau_max/N;

ind=zeros(1,length(tau));

k=zeros(1,N);

for i=0:1:N-1            %%% k is the decreasing ordered threshold's vector
    k(i+1)=tau_max-(d_tau/2+d_tau*i);
end

 for j=1:1:length(tau)     %%% for a fixed ni, scan all the threshold
    for i=1:1:N
        if tau(j)>=k(i)
            ind(j)=N+2-i;
            break
        end
    end
    if tau(j)<k(i)
        ind(j)=1;
    end
 end
end