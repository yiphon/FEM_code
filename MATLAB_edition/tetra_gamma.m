function y = tetra_gamma( gamma,node_coordinate )
%TETRA_GAMMA Summary of this function goes here
%   Detailed explanation goes here
%   gamma¡£
index=[1,2,3,4,1,2,3];
x=node_coordinate([1,4,7,10]);
z=node_coordinate([3,6,9,12]);
for i=1:1:4
    cdet=[1,x(index(i+1)),z(index(i+1));
        1,x(index(i+2)),z(index(i+2));
        1,x(index(i+3)),z(index(i+3))];
    gamma(i)=(-1)^(i+1)*det(cdet);
end
y=gamma;
end

