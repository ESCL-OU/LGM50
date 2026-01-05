function [c_ave] = calc_c_ave(c)

    c_size = size(c);
    Nr = c_size(1)-1;
    c_ave = zeros(c_size(2),1);
    
    for r = 1:Nr
        r_bar_upper = r/Nr;
        r_bar_lower = (r - 1)/Nr;
        c_ave = c_ave + (4*pi/3)*(r_bar_upper^3-r_bar_lower^3)...
            *1/2*(c(r+1, :) + c(r, :));
    end 
    c_ave = c_ave / (4/3*pi);