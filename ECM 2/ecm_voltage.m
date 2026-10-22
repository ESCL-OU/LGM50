function V = ecm_voltage(x, I, ecm_p, ocv_p, n_order)
    % z = min(max(x(:,1), 0), 1);  % clamp SOC
    z = x(1, :);

    R0 = ecm_p(1);
    OCV = polyval(ocv_p.', z);

    % Compute simulated terminal voltage
    V = OCV - R0 * I;

    % Add Vc
    for i = 1:n_order
        Vi = x(i+1, :);
        V = V + Vi;
    end

end 