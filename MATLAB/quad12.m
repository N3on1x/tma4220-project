function I = quad12(a,b,Nq,g)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
qr = {0                         , 2;
     [-sqrt(1/3), sqrt(1/3)]    , [1, 1];
     [-sqrt(3/5), 0 , sqrt(3/5)], [5/9, 8/9, 5/9];
     [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7),...
      sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)],...
     [(18-sqrt(30))/36, (18+sqrt(30))/36,...
      (18+sqrt(30))/36, (18-sqrt(30))/36]};

[t_p, w] = qr{Nq,:};

[m,n] = size(g(a));
I = zeros(m,n);

for i = 1:Nq
    I = I + g((b-a)/2 * t_p(i) + (a+b)/2)*w(i);
end
I = (b-a)/2*I;
end