function I = quadrature3D(p1,p2,p3,p4,Nq,g)
%QUADRATURE3D Evaluates surface integrals over a triangle in 2D
%   p1,p2,p2,p4: Coordinates [x,y,z] of the corners in the tetrahedron
%   Nq:       Number of barycentric quadrature points to use. Can be 1,4,5
%   g:        The integrand. A function or a matrix of functions of x and y


qr = {[1/4, 1/4, 1/4, 1/4], 1;
       0, 0;
       0, 0;
      [0.5854102, 0.1381966, 0.1381966, 0.1381966;...
       0.1381966, 0.5854102, 0.1381966, 0.1381966;...
       0.1381966, 0.1381966, 0.5854102, 0.1381966;...
       0.1381966, 0.1381966, 0.1381966, 0.5854102],...
      [0.25, 0.25, 0.25, 0.25];
      [1/4, 1/4, 1/4, 1/4; 1/2, 1/6, 1/6, 1/6;...
       1/6, 1/2, 1/6, 1/6; 1/6, 1/6, 1/2, 1/6;...
       1/6, 1/6, 1/6, 1/2],...
      [-4/5, 9/20, 9/20, 9/20, 9/20]
     };

[lb, w] = qr{Nq, :};
p = [p1' p2' p3' p4'];
[m,n] = size(g(p1(1),p1(2),p1(3)));
I = zeros(m,n);
for q = 1:Nq
    I = I + g(lb(q,:)*(p(1,:)'), lb(q,:)*(p(2,:)'), lb(q,:)*(p(3,:)'))*w(q);
end

volume = abs(det([(p2-p1)' (p3-p1)' (p4-p1)']))/6; % Volume of a tetrahedron
%volume = abs(det([p; ones(1,4)]))/6;
I = volume*I;
end

