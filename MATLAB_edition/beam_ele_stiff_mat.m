function y = beam_ele_stiff_mat( ele_info )
%BEAM_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   ���㵥Ԫ�նȾ�������ele_infoΪ[E;I;L]��ʽ��
e=ele_info(1);
i=ele_info(2);
l=ele_info(3);
kk=e*i/l^3*[12,6*l,-12,6*l;6*l,4*l^2,-6*l,2*l^2;-12,-6*l,12,-6*l;6*l,2*l^2,-6*l,4*l^2];
y=kk;
end

