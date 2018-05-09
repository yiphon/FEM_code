function y = tri_ele_b( a,node_coordinate )
%TRI_ELE_B Summary of this function goes here
%   Detailed explanation goes here
%   生成单元的B矩阵。a为面积，n_c六行一列。
x1=node_coordinate(1);
y1=node_coordinate(2);
x2=node_coordinate(3);
y2=node_coordinate(4);
x3=node_coordinate(5);
y3=node_coordinate(6);
b1=y2-y3;
b2=y3-y1;
b3=y1-y2;
c1=x3-x2;
c2=x1-x3;
c3=x2-x1;
b=1/(2*a)*[b1,0,b2,0,b3,0;0,c1,0,c2,0,c3;c1,b1,c2,b2,c3,b3];
y=b;
end

