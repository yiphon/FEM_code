function y = brick_coeffs_ui( u_expression,u )
%BRICK_COEFFS_UI Summary of this function goes here
%   Detailed explanation goes here
%   ��Ni���ǰ�u�ı��ʽ��u�еı�������ϵ������ȡϵ����ɡ�

syms u1 u2 u3 u4 u5 u6 u7 u8 s t zi

y=coeffs(u_expression,u');
y=y(end:-1:1);
y=y';

end

