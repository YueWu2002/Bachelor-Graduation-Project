function [x, y, z] = node_lobatto2D(m, plot_flag)
% [x, y, z] = node_lobatto2D(m, plot_flag)
%   compute the barycentric coordinates of the Lobatto points in the triangle
%   property: low Lebesgue number; be 1D Lobatto on edges; fully symmetric; easy generalization to higher dimensions
%
% input: 
%   m:          polynomial order
%   plot_flag:  whether to plot a figure for these points
%
% output: 
%   x, y, z:    barycentric coordinates, size = [(m+1)*(m+2)/2, 1]
%
% references:
% [1] M. G. Blyth and C. Pozrikidis. “A Lobatto interpolation grid over the triangle”. In: IMA Journal of Applied Mathematics 71.1 (Feb. 2006), pp. 153–169. issn: 0272-4960.

% checked.

N = ((m+1)*(m+2))/2;
x = zeros(N, 1);
y = zeros(N, 1);
z = zeros(N, 1);

[v, ~] = gausslobatto(m+1);
v = 0.5*(v + 1.0);

% for i = 1:m+1, j = 1: m+2-i, (k=m+3-i-j): 
% x = (1/3)*(1.0 + 2*v(i) - v(j) - v(k));
% y = (1/3)*(1.0 + 2*v(j) - v(i) - v(k));

pos = 0;
for i = 1: m+1
    for j = 1: m+2-i
        k = m+3-i-j;
        pos = pos + 1;
        x(pos) = (1/3)*(1.0 + 2*v(i) - v(j) - v(k));
        y(pos) = (1/3)*(1.0 + 2*v(j) - v(k) - v(i));
        z(pos) = (1/3)*(1.0 + 2*v(k) - v(i) - v(j));
    end
end

if (plot_flag)
    figure(1);
    x0 = [1, 0, 0];
    y0 = [0, 1, 0];
    e1=plot([x0(2),x0(3)],[y0(2),y0(3)]);
    hold on;
    e2=plot([x0(3),x0(1)],[y0(3),y0(1)]);
    hold on;
    e3=plot([x0(1),x0(2)],[y0(3),y0(2)]);
    hold on;
    plot(x0(1)*x + x0(2)*y + x0(3)*z, y0(1)*x + y0(2)*y + y0(3)*z, 'o');
    hold on;
    labels = cell(1,N);
    for k = 1: N
        labels{k} = num2str(k);
    end
    h=text(x0(1)*x + x0(2)*y + x0(3)*z, y0(1)*x + y0(2)*y + y0(3)*z, labels);
    grid off;
    axis equal;
    axis([-0.05, 1.05, -0.05, 1.05]);
    legend([e1,e2,e3],{'edge1','edge2','edge3'});
    print(gcf, '-depsc2', '-vector', 'Lobatto2D_ref', '-r300');
    hold off;

    figure(2);
    x0 = [1, 0.5, 0];
    y0 = [0, 0.5*sqrt(3.0), 0];
    plot([x0(:); x0(1)], [y0(:); y0(1)], '-');
    hold on;
    plot(x0(1)*x + x0(2)*y + x0(3)*z, y0(1)*x + y0(2)*y + y0(3)*z, 'o');
    hold on;
    labels = cell(1,N);
    for k = 1: N
        labels{k} = num2str(k);
    end
    text(x0(1)*x + x0(2)*y + x0(3)*z, y0(1)*x + y0(2)*y + y0(3)*z, labels);
    grid off;
    axis equal;
    axis([-0.05, 1.05, -0.05, sqrt(3)/2+0.05]);
    print(gcf, '-depsc2', '-vector', 'Lobatto2D_equi', '-r300');
    hold off;
end

end