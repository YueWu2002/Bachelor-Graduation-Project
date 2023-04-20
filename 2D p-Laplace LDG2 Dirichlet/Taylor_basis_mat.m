function [M, Dx, Dy] = Taylor_basis_mat(p)
% [M, Dx, Dy] = Taylor_basis_mat(p)
%   compute the matrices for Taylor basis on the right triangle (modal version)
%   should be exact in floating point
%   (we do not provide Minv, since M is exact, and cond(M) is too large)
% 
% input: 
%   p:      max total degree, p>=0
% 
% output: 
%   M:      mass matrix (accurate, exactly symmetric), size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]
%   Dx, Dy: differentiation matrices (exact), size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]

id2D = @(i,j,p) ((2*p+4-i)*(i-1))/2 + j;
N = ((p+1)*(p+2))/2;

M = nan(N, N); % exact, and symmetric
for i = 1: p+1
    for j = 1: p+2-i
        for k = 1: p+1
            for l = 1: p+2-k
                % int(K) x^(i+k-2)*y^(j+l-2)
                val = beta(i+k-1, j+l-1)/(i+j+k+l-2);
                M(id2D(i,j,p), id2D(k,l,p)) = val;
                M(id2D(k,l,p), id2D(i,j,p)) = val;
            end
        end
    end
end
% condition number is pretty big, so we try to avoid inverse M here.

Dx = sparse(N, N);
Dy = sparse(N, N);
for i = 1: p+1
    for j = 1: p+2-i
        % coef of derivative of x^(i-1)*y^(j-1), in each column
        if i > 1
            Dx(id2D(i-1,j,p), id2D(i,j,p)) = i-1;
        end
        if j > 1
            Dy(id2D(i,j-1,p), id2D(i,j,p)) = j-1;
        end
    end
end

end