function y = axisy_ele_stiff_mat( node_number,node_coordinate,b,d )
%AXISY_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元刚度矩阵。
syms r z
rmin=node_coordinate(1,node_number(3));
rmax=node_coordinate(1,node_number(1));
zmin=node_coordinate(2,node_number(3));
zmax=node_coordinate(2,node_number(1));
kk=sym(b'*d*b);
kk_intr=int(2*pi*r*kk,r,rmin,rmax);
kk_intrz=int(kk_intr,z,zmin,zmax);
y=kk_intrz;


end

