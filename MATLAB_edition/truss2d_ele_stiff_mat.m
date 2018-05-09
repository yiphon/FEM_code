function a = truss2d_ele_stiff_mat(E, A, node_coordinate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   ������Ԫ�նȾ���node_coordinateΪ��Ԫ������յ����꼯�ϡ�
x1=node_coordinate(1);
y1=node_coordinate(2);
x2=node_coordinate(3);
y2=node_coordinate(4);
l=sqrt((x2-x1)^2+(y2-y1)^2);
c=(x2-x1)/l;
s=(y2-y1)/l;
aa=E*A/l*[c^2 c*s; c*s s^2];
a=[aa -aa; -aa aa];
end

