function fun = lst_n( l )
%LST_N Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元的插值函数Ni，利用面积坐标l.
syms x y
l1=l(1);
l2=l(2);
l3=l(3);
n=sym(zeros(6,1));
n(1)=2*l1^2-l1;
n(2)=2*l2^2-l2;
n(3)=2*l3^2-l3;
n(4)=4*l1*l2;
n(5)=4*l2*l3;
n(6)=4*l3*l1;
fun=n;
end

