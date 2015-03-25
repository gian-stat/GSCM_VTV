function [ ind ] = Doppler_index( ni, N_d, ni_max )
% index compute the level's index associated to a doppler shift
d_ni=2*ni_max/N_d;

ind=zeros(1,length(ni));

k=zeros(1,N_d);

for i=0:1:N_d-1           %%% k is the decreasing ordered threshold's vector
    k(i+1)=ni_max-(d_ni/2+d_ni*i);
end

 for j=1:1:length(ni)     %%% for a fixed ni, scan all the threshold
    for i=1:1:N_d
        if ni(j)>k(i)
            ind(j)=N_d+2-i;
            break
        end
    end
    if ni(j)<k(i)
        ind(j)=1;
    end
 end
end