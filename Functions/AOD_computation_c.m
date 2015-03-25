function [ AOD, AOD_p] = AOD_computation_c( x_a, y_a, x_b, y_b, v_xt, v_yt,theta )
% Compute the AOD
% NOTE: The AOD is always calculated starting from the speed vector direction

AOD_n = atan (abs((y_b-y_a))/abs((x_b-x_a))); % (rad)

    if (y_b-y_a)>0
        if (x_b-x_a)>0
            AOD_p=AOD_n;
        else AOD_p=pi-AOD_n;
        end
    else
        if (x_b-x_a)>0
            AOD_p=2*pi-AOD_n;
        else AOD_p=pi+AOD_n;
        end
    end
    % this if is used to distinguish beetween the two directions
   if theta==0
    if v_xt>0 
        AOD=AOD_p;
    else
        if v_xt==0 % crossing case
            if v_yt>0
                 AOD=-pi/2+AOD_p;
            else AOD=pi/2+AOD_p;
            end       
        else AOD=pi+AOD_p;
        end
    end
   else
       AOD=AOD_p+(pi/2-theta);
   end
end

