function I = quadrature2D(p1,p2,p3,Nq,g)
%QUADRATURE2D Summary of this function goes here
%   Detailed explanation goes here

qr = {[1/3 1/3 1/3], 1;
      0, 0;
      [1/2 1/2 0; 1/2 0 1/2; 0 1/2 1/2], [1/3 1/3 1/3];
      [1/3, 1/3, 1/3; 3/5, 1/5, 1/5; 1/5, 3/5, 1/5; 1/5, 1/5, 3/5],...
      [-9/16, 25/48, 25/48, 25/48]
     };

[lb, w] = qr{Nq, :};
p = [p1' p2' p3'];
[m,n] = size(g(p1(1),p1(2)));
I = zeros(m,n);

for q = 1:Nq
    I = I + g(lb(q,:)*(p(1,:)'), lb(q,:)*(p(2,:)'))*w(q);
end

area = abs(det([(p2-p1)' (p3-p1)']))/2; % area of a triangle
I = area*I;
end

