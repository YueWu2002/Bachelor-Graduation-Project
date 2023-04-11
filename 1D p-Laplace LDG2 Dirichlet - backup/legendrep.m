function [y, dy] = legendrep(x, n)
% [y, dy] = legendrep(x, n)
%   compute values and derivatives of the standard Legendre polynomials
% 
% input: 
%   x:  evaluation points (can have arbitary shape)
%   n:  order of the Legendre polynomial
% 
% output:
%   y:  function values (the same shape as x)
%   dy: derivative values (the same shape as x)

% WARNING: Do NOT modify any of the code below, otherwise unexpected degeneracy in accuracy might occur!

if (n < 0)
    y = zeros(size(x));
    dy = zeros(size(x));
    return;
elseif (n == 0)
    y = ones(size(x));
    dy = zeros(size(x));
    return;
else
    % k = 0
    t1 = zeros(size(x));    % p_{-2}(x)
    t2 = zeros(size(x));    % p_{-1}(x)
    y = ones(size(x));      % p_{0}(x)
    % initialize dy
    if (mod(n, 2) ~= 0)
        dy = y;
    else
        dy = zeros(size(x));
    end

    for k = 1: n
        % compute p_{k}(x), and contribute to p_{k}'(x)
        t1 = t2;
        t2 = y;
        alpha = 1.0 / k;
        y = ((2.0-alpha)*x).*t2 - (1.0-alpha)*t1;
        % if we use "y = ((2*k-1)*x.*t2 - (k-1)*t1)/k;" the accuracy will be degenerated!
        % if we use "y = (2.0-alpha)*(x.*t2) - (1.0-alpha)*t1;" the accuracy will be degenerated!
        % if we use "y = (2*x.*t2 - t1) - alpha*(x.*t2 - t1);" the accuracy will be degenerated!
        if (mod(n-k, 2) ~= 0)
            dy = dy + (2*k+1)*y;
        end
    end
end

end