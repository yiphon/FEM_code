function y = rect_coeffs_ui( u_expression,u )
%RECT_COEFFS_UI Summary of this function goes here
%   Detailed explanation goes here
%   ��Ni���ǰ�u�ı��ʽ��u�еı�������ϵ������ȡϵ����ɡ�

syms u1 u2 u3 u4 s t

y=coeffs(u_expression,u');
y=y(end:-1:1);
y=y';

end

