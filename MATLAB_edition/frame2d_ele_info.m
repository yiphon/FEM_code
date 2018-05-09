function y = frame2d_ele_info( node_coordinate )
%FRAME2D_ELE_INFO Summary of this function goes here
%   Detailed explanation goes here
%   利用每个单元的节点坐标列阵，把单元长度、正余弦补充完整。
eles=size(node_coordinate,2);
a=zeros(3,eles);
for i=1:1:eles
    x1=node_coordinate(1,i);
    y1=node_coordinate(2,i);
    x2=node_coordinate(3,i);
    y2=node_coordinate(4,i);
    l=sqrt((x2-x1)^2+(y2-y1)^2);
    a(1,i)=l;
    c=(x2-x1)/l;
    s=(y2-y1)/l;
    a(2,i)=c;
    a(3,i)=s;
end
y=a;
end

