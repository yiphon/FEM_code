function y = frame2d_ele_stiff_mat( r,k_ele )
%FRAME2D_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   �����������ϵ�µ�Ԫ�նȾ���r��k�Ƕ�Ӧ��Ԫ��ת�ƾ���;ֲ����굥Ԫ�նȾ���
a=r'*k_ele*r;
y=a;
end

