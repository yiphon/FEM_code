function y = axisy_ele_info_for_n( node_info_for_n,node_number )
%AXISY_ELE_INFO_FOR_N Summary of this function goes here
%   Detailed explanation goes here
%   ��ÿ����Ԫ��Ӧ�Ľڵ��[1;ri;zi;ri*zi]�ļ��ϣ���������״��������Nʱ�á�
a=zeros(4,4);
for i=1:1:4
    a(i,:)=node_info_for_n(node_number(i),:);
end
y=a;
end


