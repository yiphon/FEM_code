function fun = quad2_stress_center( stress )
%QUAD2_STRESS_CENTER Summary of this function goes here
%   Detailed explanation goes here
%   生成各单元中心的应力值。
syms s t
fun=zeros(3,1);
for i=1:3
    fun(i)=double(subs(stress(i),{s,t},{0,0}));
end
end

