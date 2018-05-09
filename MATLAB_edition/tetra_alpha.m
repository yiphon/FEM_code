function y = tetra_alpha( alpha,node_coordinate )
%TETRA_ALPHA Summary of this function goes here
%   Detailed explanation goes here
%   生成单元的alpha参量。
index=[1,2,3,4,1,2,3];
x=node_coordinate([1,4,7,10]);
y=node_coordinate([2,5,8,11]);
z=node_coordinate([3,6,9,12]);
for i=1:1:4
    adet=[x(index(i+1)),y(index(i+1)),z(index(i+1));
        x(index(i+2)),y(index(i+2)),z(index(i+2));
        x(index(i+3)),y(index(i+3)),z(index(i+3))];
    alpha(i)=(-1)^(i+1)*det(adet);
end
y=alpha;
end

