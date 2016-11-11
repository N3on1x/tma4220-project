% Problem 3c

function I = line_quadrature(a,b,Nq,g)
%LINE_QUADRATURE Integrate over straight lines
%   a:  Point in R^n
%   b:  End point in R^n
%   Nq: Number of quadrature points. Can be 1,2,3,4
%   g:  Function of x, where x is a vector of n variables
qr = {0                         , 2;
     [-sqrt(1/3), sqrt(1/3)]    , [1, 1];
     [-sqrt(3/5), 0 , sqrt(3/5)], [5/9, 8/9, 5/9];
     [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7),...
      sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)],...
     [(18-sqrt(30))/36, (18+sqrt(30))/36,...
      (18+sqrt(30))/36, (18-sqrt(30))/36]};

[t_p, w] = qr{Nq,:};

r = @(t) a + (b-a).*(1+t)/2;
[m,n] = size(g(a));
I = zeros(m,n);

for i = 1:Nq
    I = I + g(r(t_p(i)))*w(i);
end
I = norm(b-a)/2*I;