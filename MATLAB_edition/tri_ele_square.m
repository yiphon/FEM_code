function y = tri_ele_square( node_coordinate )
%TRI_ELE_SQUARE Summary of this function goes here
%   Detailed explanation goes here
%   求各单元的面积。node_coordinate六行一列。
x1=node_coordinate(1);
y1=node_coordinate(2);
x2=node_coordinate(3);
y2=node_coordinate(4);
x3=node_coordinate(5);
y3=node_coordinate(6);
a=  [1,1,1;
    x1,x2,x3;
    y1,y2,y3];
b=0.5*det(a);
y=b;
end

