function y = axisy_node_info_for_n( node_coordinate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   用来求每个节点的[1;ri;zi;ri*zi]，用以求形状函数矩阵N时用，n_c=[ri;zi]。
r=node_coordinate(1);
z=node_coordinate(2);
node_info=[1,r,z,r*z];
y=node_info;
end

