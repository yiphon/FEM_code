function aa = brick_xyz_eles( n,node_coordinate )
%BRICK_XYZ_ELES Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元插值后以自然坐标表示的x\y\z变量。
syms s t zi
x=sym(zeros(8,1));
y=sym(zeros(8,1));
z=sym(zeros(8,1));
for i=1:1:8
    x(i)=node_coordinate(3*i-2);
    y(i)=node_coordinate(3*i-1);
    z(i)=node_coordinate(3*i);
end
aa=sym(zeros(3,1));
aa(1)=n'*x;
aa(2)=n'*y;
aa(3)=n'*z;
end

