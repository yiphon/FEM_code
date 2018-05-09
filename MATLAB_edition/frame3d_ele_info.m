function y = frame3d_ele_info( node_coordinate )
%FRAME3D_ELE_INFO Summary of this function goes here
%   Detailed explanation goes here
%   利用每个单元的节点坐标列阵，把单元长度、正余弦补充完整。
eles=size(node_coordinate,2);
a=zeros(4,eles);
for i=1:1:eles
    x1=node_coordinate(1,i);
    y1=node_coordinate(2,i);
    z1=node_coordinate(3,i);
    x2=node_coordinate(4,i);
    y2=node_coordinate(5,i);
    z2=node_coordinate(6,i);
    l=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    a(1,i)=l;
    cxx=(x2-x1)/l;
    cyx=(y2-y1)/l;
    czx=(z2-z1)/l;
    a(2,i)=cxx;
    a(3,i)=cyx;
    a(4,i)=czx;
end
y=a;
end

