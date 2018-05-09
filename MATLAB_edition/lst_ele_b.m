function fun = lst_ele_b( n )
%LST_ELE_B Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元B矩阵。
syms x y
ni_x=diff(n,x);
ni_y=diff(n,y);
fun=sym([ni_x,0;0,ni_y;ni_y,ni_x]);

end

