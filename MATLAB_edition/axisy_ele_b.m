function y = axisy_ele_b( jacob_eles,n )
%AXISY_ELE_B Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元的B矩阵，n一行四列。
syms r z
b=sym(zeros(4,8));
for i=1:1:4
    b(:,2*i-1:2*i)=sym([diff(n(i),r),0;n(i)/r,0;0,diff(n(i),z);diff(n(i),z),diff(n(i),r)]);
end
y=b;
end


