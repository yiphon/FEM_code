function y = axisy_ele_uv_sym( u_sym,node_number )
%AXISY_ELE_UV_SYM Summary of this function goes here
%   Detailed explanation goes here
%   ���ɸ���ԪΪ���N�к���ϵ�����õĸ��ڵ��������sym�ࡣ
a=sym(zeros(4,1));
for i=1:1:4
    a(i)=u_sym(node_number(i));
end
y=a;
end

