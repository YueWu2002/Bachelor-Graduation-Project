function [M, Dx, Dy] = Bernstein2D_basis_mat(p)
% [M, Dx, Dy] = Bernstein2D_basis_mat(p)
%   compute the matrices for Bernstein basis on the right triangle (modal version)
%   should be accurate in floating point
%   (we do not provide Minv, since M is exact, and cond(M) is too large)
% 
% input: 
%   p:      max total degree, p>=0
% 
% output: 
%   M:      mass matrix (accurate, exactly symmetric), size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]
%   Dx, Dy: differentiation matrices (exact), size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]

id2D = @(i,j,p) ((2*p+4-i)*(i-1))/2 + j; % sequence: i from 1 to p+1, j from 1 to p+2-i
N = ((p+1)*(p+2))/2;

% mass matrix
M = nan(N,N); % accurate, exactly symmetric
for i1 = 1: p+1
    for j1 = 1: p+2-i1
        k1 = p+3-i1-j1;

        for i2 = 1: p+1
            for j2 = 1: p+2-i2
                k2 = p+3-i2-j2;

                val = Bernstein2D_inner_product(i1-1,j1-1,k1-1, i2-1,j2-1,k2-1);
                M(id2D(i1,j1,p), id2D(i2,j2,p)) = val;
            end
        end
    end
end
% condition number is pretty big, so we try to avoid inverse M here.

% differential relations
Dx = sparse(N,N);
Dy = sparse(N,N);
for i = 1: p+1
    for j = 1: p+2-i
        k = p+3-i-j;

        % coef of derivatives of P_{i-1,j-1,k-1}, in each column
        Dx(id2D(i,j,p), id2D(i,j,p)) = (i-1) - (k-1);
        Dy(id2D(i,j,p), id2D(i,j,p)) = (j-1) - (k-1);

        if i-1 > 0
            Dx(id2D(i-1,j+1,p), id2D(i,j,p)) = (j-1) + 1;
            Dx(id2D(i-1,j,p), id2D(i,j,p)) =  (k-1) + 1;
        end

        if j-1 > 0
            Dy(id2D(i+1,j-1,p), id2D(i,j,p)) = (i-1) + 1;
            Dy(id2D(i,j-1,p), id2D(i,j,p)) = (k-1) + 1;
        end
        
        if k-1 > 0
            Dx(id2D(i+1,j,p), id2D(i,j,p)) = -((i-1) + 1);
            Dx(id2D(i,j+1,p), id2D(i,j,p)) = -((j-1) + 1);

            Dy(id2D(i+1,j,p), id2D(i,j,p)) = -((i-1) + 1);
            Dy(id2D(i,j+1,p), id2D(i,j,p)) = -((j-1) + 1);
        end

    end
end


end


function [val] = Bernstein2D_inner_product(i1,j1,k1, i2,j2,k2)
% compute the inner product of Bernstein polynomials of order (i1,j1,k1) and (i2,j2,k2)
%   domain: the right triangle, with vertices (1,0), (0,1), (0,0)
% 
% input: 
%   i1,j1,k1, i2,j2,k2:     integers, scalar, >=0
% 
% output: 
%   val:    double, scalar

if min([i1,j1,k1, i2,j2,k2]) < 0
    error('Bernstein2D_inner_product: error input!');
end

n1 = i1 + j1 + k1;
n2 = i2 + j2 + k2;
n = n1 + n2;
i = i1 + i2;
j = j1 + j2;
k = k1 + k2;

% val = gamma(n1+1) * gamma(n2+1) * gamma(i+1) * gamma(j+1) * gamma(k+1) / (gamma(i1+1) * gamma(i2+1) * gamma(j1+1) * gamma(j2+1) * gamma(k1+1) * gamma(k2+1) * gamma(n1+n2+3));
% val = beta(n1+1, n2+1) / ((n+2) * ((i+1) * beta(i1+1, i2+1)) * ((j+1) * beta(j1+1, j2+1)) * ((k+1) * beta(k1+1, k2+1)));
val = double( nchoosek(int64(i), int64(i1)) * nchoosek(int64(j), int64(j1)) * nchoosek(int64(k), int64(k1)) ) / double( int64(n+2) * int64(n+1) * nchoosek(int64(n), int64(n1)) ); % avoid overflow, and increse accuracy (which proves effective in very high-order computation)

end