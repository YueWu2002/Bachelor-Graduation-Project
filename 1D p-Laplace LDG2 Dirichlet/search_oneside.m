function [x, y] = search_oneside(fun, y0, x_guess, delta)

% we re-implemeent the brute-force method
% convergence is faster, however, accuracy is not as good

if x_guess <= 0.0
    x = 0.0;
    y = y0;
    return;
end

gold = sqrt(1.25) - 0.5;
list = gold.^linspace(1.0, -1.0, 3);

y_trial = zeros(length(list), 1);
x = x_guess;

while true
    for k = 1: length(list)
        y_trial(k) = fun(list(k)*x);
    end
    [y_min, id] = min(y_trial(:));
    if y_min >= y0 && abs(x) > delta
        x = 0.25*x;
    else
        x = list(id(1))*x;
        y = y_min;
        return;
    end
end


end