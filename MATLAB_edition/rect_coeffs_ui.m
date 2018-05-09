function y = rect_coeffs_ui( u_expression,u )
%RECT_COEFFS_UI Summary of this function goes here
%   Detailed explanation goes here
%   求Ni，是把u的表达式以u中的变量整理系数，提取系数组成。

syms u1 u2 u3 u4 s t

y=coeffs(u_expression,u');
y=y(end:-1:1);
y=y';

end

