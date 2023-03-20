function [y1] = interpolate(x0, y0, x1)
% size of x0 and y0 must agree
% y1 will have the same shape as x1

ss = size(x1);
x1 = reshape(x1, [prod(ss), 1]);
y1 = zeros(prod(ss), 1);

nn = numel(x0);
x0 = reshape(x0, [nn, 1]);
y0 = reshape(y0, [nn, 1]);

for i = 2: nn
    for j = nn: -1: i
        y0(j) = (y0(j) - y0(j-1))/(x0(j) - x0(j-i+1));
    end
end

y1(:) = y0(nn);
for i = nn-1: -1: 1
    y1(:) = y1(:).*(x1(:) - x0(i)) + y0(i);
end

y1 = reshape(y1, ss);

end