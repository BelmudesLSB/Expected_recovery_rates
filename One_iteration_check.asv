%% This script performs a one iteration check on the solutions obtained from the .mex file.

Q1 = zeros(size(Q));
B1_policy = zeros(size(B_policy)); 
D1_policy = zeros(size(D_policy));
V1_r = zeros(size(V_r));
V1_d = zeros(size(V_d));
V1 = zeros(size(V));

% Update the value function:
for y = 1:params.y_grid_size
    for b = 1:params.b_grid_size
        if V_r(y,b)>= V_d(y,b)
            V1(y,b)= V_r(y,b);
            D1_policy(y,b) = 0;
        else
            V1(y,b)= V_d(y,b);
            D1_policy(y,b) = 1;
        end
    end
end

% Update price given default policy:
for y = 1:params.y_grid_size
    for b_prime = 1:params.b_grid_size
        for y_prime = 1:params.y_grid_size
            Q1(y,b_prime) =  Q1(y,b_prime) + P(y, y_prime) * (1-D1_policy(y_prime, b_prime)) * (1/(1+params.r));
        end
    end
end

% Update the value at default:
for y = 1:params.y_grid_size
    for b = 1:params.b_grid_size
        E_vd = 0;
        E_v = 0;
        for y_prime = 1:params.y_grid_size
            E_v = E_v + P(y,y_prime) * V(y_prime, params.b_grid_size);
            E_vd = E_vd + P(y,y_prime) * V_d(y_prime, params.b_grid_size);
        end
    V1_d(y, b) = (1/(1-params.gamma)) * (Y_grid_default(y) + params.alpha * B_grid(b))^((1-params.gamma)) + params.beta * (params.theta * E_v + (1-params.theta) * E_vd);
    end
end

% Update value of repayment and bond policy:
for y = 1:params.y_grid_size
    for b = 1:params.b_grid_size 
        temp = -1000000000;
        for b_prime = 1:params.b_grid_size
            E_v_bprime = 0;
            for y_prime = 1:params.y_grid_size
                E_v_bprime = E_v_bprime + P(y, y_prime) * V1(y_prime, b_prime);
            end
            if Y_grid(y) + B_grid(b) - Q1(y, b_prime) * B_grid(b_prime) > 0
                util = (1/(1-params.gamma)) * (Y_grid(y) + B_grid(b) - Q1(y, b_prime) * B_grid(b_prime))^(1-params.gamma) + params.beta * E_v_bprime;
                if util > temp
                    temp = util;
                    B1_policy(y,b) = b_prime;
                end
            end
        end
        V1_r(y,b) = temp;
    end
end


max(max(abs(V1 - V)))
max(max(abs(V1_d - V_d)))
max(max(abs(V1_r - V_r)))
max(max(abs(Q1 - Q)))
max(max(abs(D1_policy-D_policy)))
max(max(abs(B_policy(D1_policy<1) - B1_policy(D1_policy<1))))


    
