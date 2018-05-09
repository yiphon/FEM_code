function y = tetra_delta(delta,node_coordinate)
%TETRA_DELTA Summary of this function goes here
%   Detailed explanation goes here
%   delta.
index=[1,2,3,4,1,2,3];
x=node_coordinate([1,4,7,10]);
y=node_coordinate([2,5,8,11]);
for i=1:1:4
    ddet=[1,x(index(i+1)),y(index(i+1));
        1,x(index(i+2)),y(index(i+2));
        1,x(index(i+3)),y(index(i+3))];
    delta(i)=(-1)^(i)*det(ddet);
end
y=delta;
end

