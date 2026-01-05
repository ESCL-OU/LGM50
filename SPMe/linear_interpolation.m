function vq = linear_interpolation(x, v, xq)
    % 1D Linear Interpolation
    % x  : Vector of known x values (must be sorted in ascending order)
    % v  : Vector of known y values corresponding to x
    % xq : Vector of query points for interpolation
    % vq : Interpolated values at the query points

    % Ensure x and v are column vectors
    x = x(:);
    v = v(:);
    
    % Preallocate output array
    if isa(v,'double')
        % vq = zeros(size(xq));
        vq = [];
        v_type = 'd'; % Double
    else
        % vq = casadi.MX.sym('vq', length(xq));
        vq = vertcat();
        v_type = 'c'; % CasADi
    end
    % Loop through each query point
    for i = 1:length(xq)
        % Handle extrapolation on the left
        if xq(i) <= x(1)
            % vq(i) = v(1);
            % vq(end+1) = v(1);
            vqk = v(1);
        % Handle extrapolation on the right
        elseif xq(i) >= x(end)
            % vq(i) = v(end);
            % vq(end+1) = v(end);
            vqk = v(end);
        else
            % Find the interval [x_k, x_k+1] where x_k <= xq < x_k+1
            k = find(x <= xq(i), 1, 'last');
            
            % Calculate the slope (difference quotient)
            slope = (v(k+1) - v(k)) / (x(k+1) - x(k));
            
            % Apply the linear interpolation formula
            % vq(i) = v(k) + slope * (xq(i) - x(k));
            % vq(end+1) = v(k) + slope * (xq(i) - x(k));
            vqk = v(k) + slope * (xq(i) - x(k));
        end
        if v_type == 'd'
            vq = [vq, vqk];
        else
            vq = vertcat(vq, vqk);
        end
    end
end
