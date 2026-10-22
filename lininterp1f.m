function yinterp = lininterp1f(x, y, xi, ydefault)
%LININTERP1F Lightweight MATLAB fallback for legacy SPMeT calls.
% Emulates linear interpolation with optional default fill outside bounds.

    if nargin < 4 || isempty(ydefault)
        yinterp = interp1(x, y, xi, 'linear', 'extrap');
    else
        yinterp = interp1(x, y, xi, 'linear', ydefault);
    end
end
