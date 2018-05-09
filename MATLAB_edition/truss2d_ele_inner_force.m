function a = truss2d_ele_inner_force( E, A, u, node_coordinate)
%ELE_INNER_FORCE Summary of this function goes here
%   Detailed explanation goes here
%   ���ɸ���Ԫ������uΪ��Ӧ��Ԫ���ڵ��λ�ơ�
x1=node_coordinate(1);
y1=node_coordinate(2);
x2=node_coordinate(3);
y2=node_coordinate(4);
l=sqrt((x2-x1)^2+(y2-y1)^2);
c=(x2-x1)/l;
s=(y2-y1)/l;
a=E*A/l*[-c,-s,c,s]*u;
end

