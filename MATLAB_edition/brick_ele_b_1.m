function y = brick_ele_b_1( Ni,x,y,z )
%BRICK_ELE_B_1 Summary of this function goes here
%   Detailed explanation goes here
%   生成单元的B矩阵。Ni，x，y均为s\t的函数。
syms s t zi 
b=sym(zeros(6,3));
j=det(jacobian([x,y,z],[s,t,zi]));
ni_x=det(jacobian([Ni,y,z],[s,t,zi]))/j;
ni_y=det(jacobian([x,Ni,z],[s,t,zi]))/j;
ni_z=det(jacobian([x,y,Ni],[s,t,zi]))/j;
b=[ni_x,0,0;0,ni_y,0;0,0,ni_z;ni_y,ni_x,0;0,ni_z,ni_y;ni_z,0,ni_x];
y=b;



end

