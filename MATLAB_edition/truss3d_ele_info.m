function a = truss3d_ele_info(node_coordinate)
%ELE_ANGLE_3DTRUSS Summary of this function goes here
%   Detailed explanation goes here
%   生成单元基本信息列阵。
x1=node_coordinate(1);
y1=node_coordinate(2);
z1=node_coordinate(3);
x2=node_coordinate(4);
y2=node_coordinate(5);
z2=node_coordinate(6);
b=zeros(4,1);
l=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
b(1)=(x2-x1)/l;
b(2)=(y2-y1)/l;
b(3)=(z2-z1)/l;
b(4)=l;
a=b;
end

