function I = quadr1D(a, b, Nq, g)
%QUADR1D Summary of this function goes here
%   Detailed explanation goes here

qr = {0                         , 2;
     [-sqrt(1/3), sqrt(1/3)]    , [1, 1];
     [-sqrt(3/5), 0 , sqrt(3/5)], [5/9, 8/9, 5/9];
     [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7),...
      sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)],...
     [(18-sqrt(30))/36, (18+sqrt(30))/36,...
      (18+sqrt(30))/36, (18-sqrt(30))/36]};

[t_p, w] = qr{Nq,:};

r = @(t) a + (b-a).*(1+t)/2;
I = 0;

for i = 1:Nq
    val = r(t_p(i));
    I = I + g(val(1),val(2))*w(i);
end
I = norm(b-a)/2*I;
end

