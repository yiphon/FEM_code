function y = tetra_main_stress( sigma )
%TETRA_MAIN_STRESS Summary of this function goes here
%   Detailed explanation goes here
%   生成单元的主应力与主方向。
y={};
sx=sigma(1);
sy=sigma(2);
sz=sigma(3);
txy=sigma(6);
txz=sigma(5);
tyz=sigma(4);
sigma_mat=[sx,txy,txz;
           txy,sy,tyz;
           txz,tyz,sz];
[a,b]=eig(sigma_mat);
b=[b(1,1);b(2,2);b(3,3)];
y{1,1}=b;
y{1,2}=a;
end

