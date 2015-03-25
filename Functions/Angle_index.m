function [ ind ] = Angle_index( angle, N_angle)
% index compute the level's index associated to a tau
N_angle_p=N_angle-1;

d_ang=2*pi/N_angle_p;
ind=zeros(1,length(angle));

k=zeros(1,N_angle_p);

for i=0:1:N_angle_p-1            %%% k is the decreasing ordered threshold's vector
    k(i+1)=2*pi-(d_ang/2+d_ang*i);
end
angle_p=zeros(1,length(angle));
 for j=1:1:length(angle)     %%% for a fixed angle, scan all the threshold
     if angle(j)>=2*pi
         angle_p(j)=angle(j)-2*pi;
     else
         angle_p(j)=angle(j);
     end
    for i=1:1:N_angle_p
        if angle_p(j)>=k(i)
            ind(j)=N_angle_p+2-i;
            break
        end
    end
    if angle_p(j)<k(i)
        ind(j)=1;
    end
 end
end