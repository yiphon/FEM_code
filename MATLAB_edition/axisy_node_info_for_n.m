function y = axisy_node_info_for_n( node_coordinate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   ������ÿ���ڵ��[1;ri;zi;ri*zi]����������״��������Nʱ�ã�n_c=[ri;zi]��
r=node_coordinate(1);
z=node_coordinate(2);
node_info=[1,r,z,r*z];
y=node_info;
end

