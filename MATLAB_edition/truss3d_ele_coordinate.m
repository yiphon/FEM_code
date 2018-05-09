function y = truss3d_ele_coordinate( node_coordinate,node_number )
%TRUSS3D_ELE_COORDINATE Summary of this function goes here
%   Detailed explanation goes here
%   ����ÿ����Ԫ��Ӧ�Ľڵ���������
eles=size(node_number,2);
k=zeros(6,eles);
for i=1:1:eles
    for j=1:1:2
    k(3*j-2:3*j,i)=node_coordinate(:,node_number(j,i));
    end
end
y=k;


end

