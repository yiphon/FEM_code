function y = axisy_ele_b_1( Ni,r,z )
%AXISY_ELE_B_1 Summary of this function goes here
%   Detailed explanation goes here
%   生成单元的B矩阵。Ni，r，z均为s\t的函数。
syms s t 
b=sym(zeros(4,2));
j=det(jacobian([r,z],[s,t]));
ni_r=det(jacobian([Ni,z],[s,t]))/j;
ni_z=det(jacobian([r,Ni],[s,t]))/j;
b=[ni_r,0;Ni/r,0;0,ni_z;ni_z,ni_r];
y=b;

end

