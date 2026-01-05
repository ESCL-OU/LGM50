function x_dot = ode_ecm(x, I, p, ecm_p, n_order)
    Q = p.OneC * 3600;  % Capacity in Coulombs

    R1 = ecm_p(2);
    C1 = ecm_p(3);
    x_dot_z =  -I / Q;                          % SOC dynamics
    
    x_dot = [x_dot_z];
    % Second order or higher
    for i = 1:n_order
        Ri = ecm_p(i*2);% 4 6 8
        Ci = ecm_p(i*2+1);% 5 7 9
        x_dot_i = -(1/Ci) * I - (1/(Ri*Ci)) * x(i+1);
        x_dot = [x_dot; x_dot_i];
    end
end
