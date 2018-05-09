function y = frame3d_r_mat( ele_info )
%FRAME3D_R_MAT Summary of this function goes here
%   Detailed explanation goes here
%   ���ɸ���Ԫ��ת�ƾ��󡣣�ȡ�ֲ�����ϵ��ԭ����y����XOYƽ���ڣ���ele_info[c;c;c]��
%   %һ�־ֲ�����ϵ��ȡ��ԭ���Ǿֲ�x��ȫ��Z��˵þֲ�y�����᷽�򣬼����������ȡ����
cxx=ele_info(1);
cyx=ele_info(2);
czx=ele_info(3);
if cxx==0 & cyx==0
    if czx>0
        r=[0,0,1;0,1,0;-1,0,0];
    else
        r=[0,0,-1;0,1,0;1,0,0];
    end
else
    d=sqrt(cxx^2+cyx^2);
    cxy=-cyx/d;
    cyy=cxx/d;
    czy=0;
    cxz=-czx*cyy;
    cyz=-czx*cyx/d;
    czz=d;
    r=[cxx,cyx,czx;cxy,cyy,czy;cxz,cyz,czz];
end
r0=zeros(3,3);
rr=[r,r0,r0,r0;r0,r,r0,r0;r0,r0,r,r0;r0,r0,r0,r];
y=rr;
end

