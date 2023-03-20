p_equ = 10;
q_equ = p_equ/(p_equ-1);

uu = @(x) (x + 1.0).^((p_equ-2)/(p_equ-1)); % u(x)
dduu = @(x) ((p_equ-2)/(p_equ-1)) * (x + 1.0).^(-1.0/(p_equ - 1));    % u'(x)
qq = @(x) (abs((p_equ-2)/(p_equ-1))^(p_equ-1))*sign(p_equ-2)./(x+1);    % flux[u](x)
ff = @(x) (abs((p_equ-2)/(p_equ-1))^(p_equ-1))*sign(p_equ-2)./((x+1).^2);  % f(x)

a = 0;
b = 3*pi;
gl = uu(a);
gr = uu(b);
type = 'lobatto';

p = 2;
switch type
case 'legendre'
    [qx, qw] = gausslegendre(p+1);
case 'lobatto'
    [qx, qw] = gausslobatto(p+1);
end
[~, ~, D, ~, ~, ~, ~] = matrix_nodal1D(qx);

err = [];
h = [];
N = [];
for k = 1: 6
    N(k) = 2^(k+3);
    xx = linspace(a, b, N(k)+1);
    h(k) = max(abs(xx(2:end) - xx(1:end-1)));

    [q, u, du] = pLaplace1D_LDG(ff, xx, gl, gr, p, p_equ, 1, 0.5, type, true);
    
    err(k,:) = [0, 0, 0, 0];
    for j = 1: length(xx) - 1
        xc = 0.5*(xx(j) + xx(j+1));
        for l = 1: p+1
            % error at the i-th node

            xp = xc + 0.5*(xx(j+1)-xx(j))*qx(l);
            % Lp error
            err(k,1) = err(k,1) + (xx(j+1)-xx(j))/2 * qw(l) * abs(u((j-1)*(p+1) + l) - uu(xp))^p_equ;

            % Lq error
            err(k,2) = err(k,2) + (xx(j+1)-xx(j))/2 * qw(l) * abs(q((j-1)*(p+1) + l) - qq(xp))^q_equ;

            % L2 error
            err(k,3) = err(k,1) + (xx(j+1)-xx(j))/2 * qw(l) * abs(u((j-1)*(p+1) + l) - uu(xp))^2;
            err(k,4) = err(k,2) + (xx(j+1)-xx(j))/2 * qw(l) * abs(q((j-1)*(p+1) + l) - qq(xp))^2;
        end
    end
    err(k,:) = err(k,:).^[1/p_equ, 1/q_equ, 1/2, 1/2];
end

ord = nan(size(err));
for k = 2: size(ord, 1)
    ord(k, :) = log(err(k,:)./err(k-1,:))/log(h(k)./h(k-1));
end

format short e
disp(err)

format short
disp(ord)


for i = 1: size(err, 1)
    fprintf('%d ', N(i));
    for j = 1: size(err, 2)
        fprintf('& %.4e & %.4f ', [err(i,j), ord(i,j)]);
    end
    fprintf('\n');
end