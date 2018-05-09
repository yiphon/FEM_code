function y = tetra_beta( beta,node_coordinate )
%TETRA_BETA Summary of this function goes here
%   Detailed explanation goes here
%   生成beta参量。
index=[1,2,3,4,1,2,3];
y=node_coordinate([2,5,8,11]);
z=node_coordinate([3,6,9,12]);
for i=1:1:4
    bdet=[1,y(index(i+1)),z(index(i+1));
        1,y(index(i+2)),z(index(i+2));
        1,y(index(i+3)),z(index(i+3))];
    beta(i)=(-1)^(i)*det(bdet);
end
y=beta;
end

