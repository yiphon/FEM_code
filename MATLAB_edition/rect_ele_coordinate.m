function y = rect_ele_coordinate( node_coordinate,node_number )
%RECT_ELE_COORDINATE Summary of this function goes here
%   Detailed explanation goes here
%   ����ÿ����Ԫ��Ӧ�Ľڵ���������
eles=size(node_number,2);
k=zeros(8,eles);
for i=1:1:eles
    for j=1:1:4
    k(2*j-1:2*j,i)=node_coordinate(:,node_number(j,i));
    end
end
y=k;

end

