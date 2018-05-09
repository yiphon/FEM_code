function mat = brick_ele_b( ni,x,y,z )
%BRICK_ELE_B Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元B矩阵。j_e=dni/d(s,t,zi).
syms s t zi

j=det(jacobian([x,y,z]));
ni_x_nume=det(jacobian([ni,y,z]));
ni_y_nume=det(jacobian([x,ni,z]));
ni_z_nume=det(jacobian([x,y,ni]));
ni_x=ni_x_nume/j;
ni_y=ni_y_nume/j;
ni_z=ni_z_nume/j;

mat=[ni_x,0,0;0,ni_y,0;0,0,ni_z;ni_y,ni_x,0;0,ni_z,ni_y;ni_z,0,ni_x];
end






