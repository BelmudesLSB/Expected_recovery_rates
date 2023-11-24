function [stats, simulations] = Run_Simulations(p, p_sim, Eq, Rand_Vec)
    
    % Size of the simulation:
    %TBurn = 12*10^3;
    %T = TBurn + 12*10^6;
    
    T = p_sim.T;
    TBurn = p_sim.TBurn;

    % Empty container:
    yt = NaN*zeros(T,1);
    yt_default = NaN*zeros(T,1);
    bt_lowr = NaN*zeros(T,1);
    bt_highr = NaN*zeros(T,1);
    Dt = NaN*zeros(T+1,1);

    % Generate the income path:
    current_y = round(p.y_grid_size/2);
    current_b_lowr = p.b_grid_size_lowr;
    current_b_highr = p.b_grid_size_highr;
    Dt(1) = 0;

    for step = 1:T

        yt(step) = current_y;
        next_y = randsample(p.y_grid_size, 1, true, Eq.P(current_y, :));

        if Dt(step) == 0
            bt_lowr(step) = current_b_lowr;
            bt_highr(step) = current_b_highr;        
            next_bt_lowr = Eq.B_policy_lowr(current_b_highr, current_b_lowr, current_y);
            next_bt_highr = Eq.B_policy_highr(current_b_highr, current_b_lowr, current_y);
            Dt(step+1) = Eq.D_policy(next_bt_highr, next_bt_lowr, next_y);
            current_b_lowr = next_bt_lowr;
            current_b_highr = next_bt_highr; 
            current_y = next_y;
        else
            if Rand_Vec.theta(step) > p.theta % Fail to reenter
                bt_lowr(step) = current_b_lowr;
                bt_highr(step) = current_b_highr;
                current_b_lowr = p.b_grid_size_lowr;
                current_b_highr = p.b_grid_size_highr;
                current_y = next_y;
                Dt(step+1) = 1;
            else
                bt_lowr(step) = current_b_lowr;
                bt_highr(step) = current_b_highr;
                Dt(step+1) = 0;
                current_b_lowr = p.b_grid_size_lowr;
                current_b_highr = p.b_grid_size_highr;
                current_y = next_y;
            end
        end
    end

    simulations.Default_policy = Dt(TBurn:end-1);
    simulations.B_low = Eq.B_grid_lowr(bt_lowr(TBurn:end));
    simulations.B_low_index = bt_lowr(TBurn:end);
    simulations.B_high = Eq.B_grid_highr(bt_highr(TBurn:end));
    simulations.B_high_index = bt_highr(TBurn:end);
    simulations.Y = Eq.Y_grid(yt(p_sim.TBurn:end));
    simulations.Y_index = yt(p_sim.TBurn:end);

    stats.Y = mean(simulations.Y(simulations.Default_policy == 0));
    stats.B_lowr = mean(simulations.B_low(simulations.Default_policy == 0));
    stats.B_highr = mean(simulations.B_high(simulations.Default_policy == 0));
    stats.Default_policy = sum(simulations.Default_policy)/T;

end