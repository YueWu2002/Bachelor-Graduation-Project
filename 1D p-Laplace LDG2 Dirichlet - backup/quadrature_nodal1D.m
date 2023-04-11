function [qx, qw] = quadrature_nodal1D(n, type)

switch type
case 'legendre'
    [qx, qw] = gausslegendre(n);
case 'lobatto'
    [qx, qw] = gausslobatto(n);
case 'equid'
    [qx] = linspace(-1, 1, n)';
    [qw] = quadweight_nodal1D(qx);
otherwise
    disp('Wrong collocation type!');
end

end

function [w] = quadweight_nodal1D(x)
% directly solve the weights using Taylor basis
n = length(x);
b = zeros(n, 1);
b(1:2:n) = (2.0)./(1:2:n);
w = fliplr(vander(x))'\b;
% remark: since the quadrature points have error, the integration accuracy will be the best if we use the above method to compute the weights, given n is not too big.
end