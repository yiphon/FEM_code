function y = axisy_n_coeffs_for_rz( ele_u_sym,ele_info_for_n )
%AXISY_N_COEFFS_FOR_RZ Summary of this function goes here
%   Detailed explanation goes here
%   ���ɸ���Ԫ��������λ�Ʋ�ֵ��ʽ��r\zΪ����õ���ϵ���������Ǹ�˹��Ԫ����
a=sym(ele_info_for_n);
b=ele_u_sym;
c=a\b;
y=c;


end

