function y = axisy_n_coeffs_for_rz( ele_u_sym,ele_info_for_n )
%AXISY_N_COEFFS_FOR_RZ Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元求解出来的位移插值形式以r\z为整理得到的系数，本质是高斯消元法。
a=sym(ele_info_for_n);
b=ele_u_sym;
c=a\b;
y=c;


end

