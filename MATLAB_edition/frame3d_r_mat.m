function y = frame3d_r_mat( ele_info )
%FRAME3D_R_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元的转移矩阵。（取局部坐标系的原则是y轴在XOY平面内？）ele_info[c;c;c]。
%   %一种局部坐标系的取法原则是局部x与全局Z叉乘得局部y坐标轴方向，即本函数里的取法。
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

