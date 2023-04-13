function [x, y] = golden_search_oneside(fun, y0, x_guess, delta)
% one-sided golden section search, which can minimize the number of 'fun' call
% 
% input:
%   fun:        convex function on [0, +inf), with f(+inf) = +inf
%   y0:         fun evaluated at 0
%   x_guess:    initial guess, x_guess > 0
%   delta:      accuracy control for x, should be assigned relatively with x_guess
% 
% output: 
%   x:          minimal point, scalar, x >= 0
%   y:          minimal value, scalar, y <= y0
% 
% REMARK: it gives x = x_guess if x_guess is the exact minimal point (consistency)

gold = sqrt(1.25) - 0.5;

if x_guess <= 0.0
    x = 0.0;
    y = y0;
    return;
end

y_at_guess = fun(x_guess);

if y0 <= y_at_guess
    % search in the interval [0, x_guess] with boundary values [y0, y_at_guess]

    % make them into arrays with length 3
    xx = [0.0, (1.0-gold)*x_guess, x_guess];
    yy = [y0, fun(xx(2)), y_at_guess];

    % accelerate convergence (we keep that y0 <= yy(3))
    % maintain length 3, and maintain the golden section structure of xx
    while xx(end) > delta
        if y0 <= yy(2)
            % should further shrink the interval
            xx = [0.0, (1.0-gold)*xx(2), xx(2)];
            yy = [y0, fun(xx(2)), yy(2)];
        else % yy(2) < y0 <= yy(3) 
            break;
        end
    end

    xx = [0.0, xx(2), gold*xx(2) + (1.0-gold)*xx(3), xx(3)];
    yy = [y0, yy(2), fun(xx(3)), yy(3)];

else % y0 > y_at_guess

    % make them into arrays with length 3
    xx = [0.0, x_guess, x_guess/gold];
    yy = [y0, y_at_guess, fun(xx(3))];
    
    % maintain length 3, and maintain the golden section structure of xx
    while yy(2) > yy(3)
        xx = [xx(2:3), xx(3)/gold];
        yy = [yy(2:3), fun(xx(3))];
    end
    % when the iteration stops: yy(1) > yy(2) <= yy(3)
    % search in the interval [xx(1), xx(3)] with boundary values [yy(1), yy(3)]

    xx = [xx(1), xx(2), gold*xx(2) + (1.0-gold)*xx(3), xx(3)];
    yy = [yy(1), yy(2), fun(xx(3)), yy(3)];

end


% in-interval golden section search for convex functions
% the minimal point is guaranteed to be included in the interval

while xx(end) - xx(1) > delta
    if yy(2) > yy(3) % must be strict to make x close to 0
        xx = [xx(2), xx(3), gold*xx(3) + (1.0-gold)*xx(4), xx(4)];
        yy = [yy(2), yy(3), fun(xx(3)), yy(4)];
    else
        xx = [xx(1), (1.0-gold)*xx(1) + gold*xx(2), xx(2), xx(3)];
        yy = [yy(1), fun(xx(2)), yy(2), yy(3)];
    end
end

% search from left to right
id = 1;
minval = yy(1);
for k = 2: length(yy)
    if yy(k) < minval % must be strict to accelerate convergence
        id = k;
        minval = yy(k);
    end
end
x = xx(id);
y = minval;

end