function y = tri_ele_coordinate( node_coordinate,node_number )
%TRI_ELE_COORDINATE Summary of this function goes here
%   Detailed explanation goes here
%   ����ÿ����Ԫ��Ӧ�Ľڵ���������
eles=size(node_number,2);
k=zeros(6,eles);
for i=1:1:eles
    for j=1:1:3
    k(2*j-1:2*j,i)=node_coordinate(:,node_number(j,i));
    end
end
y=k;

end

