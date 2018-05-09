function y = tetra_ele_b( v,beta,gamma,delta )
%TETRA_ELE_B Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元B矩阵。
b=1/(6*v)*[beta,0,0;0,gamma,0;0,0,delta;
    gamma,beta,0;0,delta,gamma;delta,0,beta];
y=b;
end

