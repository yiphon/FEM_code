function y = axisy_ele_uv_sym( u_sym,node_number )
%AXISY_ELE_UV_SYM Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元为求解N中函数系数设置的各节点坐标矩阵，sym类。
a=sym(zeros(4,1));
for i=1:1:4
    a(i)=u_sym(node_number(i));
end
y=a;
end

