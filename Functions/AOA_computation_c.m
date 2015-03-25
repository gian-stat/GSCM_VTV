function [ AOA, AOA_p ] = AOA_computation_c( x_a, y_a, x_b, y_b, v_xr, v_yr,theta )
% Compute the AOA
% NOTE: The AOA is always calculated starting from the speed vector direction

AOA_n = atan (abs((y_b-y_a))/abs((x_b-x_a))); % (rad)

    if (y_b-y_a)>0
        if (x_b-x_a)>0
            AOA_p=pi+AOA_n;
        else AOA_p=2*pi-AOA_n;
        end
    else
        if (x_b-x_a)>0
            AOA_p=pi-AOA_n;
        else AOA_p=AOA_n;
        end
    end
    
    % this if is used to distinguish beetween the two directions
    if theta==0
    if v_xr>0 
        AOA=AOA_p;
    else
        if v_xr==0 % crossing case
            if v_yr>0
                 AOA=-pi/2+AOA_p;
            else AOA=pi/2+AOA_p;
            end       
        else AOA=pi+AOA_p;
        end
    end
    else
       AOA=AOA_p+(pi/2-theta);
   end
end

