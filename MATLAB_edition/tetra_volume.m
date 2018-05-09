function y = tetra_volume( node_coordinate )
%TETRA_VOLUME Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元的体积。
x=node_coordinate([1,4,7,10]);
y=node_coordinate([2,5,8,11]);
z=node_coordinate([3,6,9,12]);
vdet=[1,x(1),y(1),z(1);
    1,x(2),y(2),z(2);
    1,x(3),y(3),z(3);
    1,x(4),y(4),z(4)];
v=1/6*det(vdet);
y=v;
end

