function a = truss3d_ele_inner_force( E,A,ele_info,u)
%ELE_INNER_FORCE_3DTRUSS Summary of this function goes here
%   Detailed explanation goes here
%   求单元桁架中的内力。
cx=ele_info(1);
cy=ele_info(2);
cz=ele_info(3);
l=ele_info(4);
a=E*A/l*[-cx,-cy,-cz,cx,cy,cz]*u;

end

