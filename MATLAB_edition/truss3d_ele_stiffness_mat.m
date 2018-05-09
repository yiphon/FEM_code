function a =truss3d_ele_stiffness_mat( E,A,ele_info )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   生成单元刚度矩阵，利用单元信息矩阵。
cx=ele_info(1);
cy=ele_info(2);
cz=ele_info(3);
l=ele_info(4);
aa=E*A/l*[cx^2,cx*cy, cx*cz;cy*cx,cy^2,cy*cz;cz*cx,cz*cy,cz^2];
a=[aa,-aa;-aa,aa];
end

