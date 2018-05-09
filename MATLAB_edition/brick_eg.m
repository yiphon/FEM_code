%六面体元算例。
%也可以用来构建等参元。
source=xlsread('brick_eg_input.xlsx');
source(all(isnan(source),2),:)=[];
eles=source(1,1);
nodes=source(1,2);
e=source(2,1);
nu=source(2,2);
node_number=source(3:2+eles,1:8)';   %每个单元对应节点的编号。
node_coordinate=source(3+eles:2+eles+nodes,1:3)';
node_coordinate=brick_ele_coordinate(node_coordinate,node_number);   %每个单元对应节点的坐标。
unknown_u_index=source(3+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';       %unknown_u_index是U列向量中未知位移位置的索引。
known_f=source(4+eles+nodes,1:size(unknown_u_index,1))';     %f_known为已知的节点力，包含等效节点载荷。
known_u=source(5+eles+nodes,1:3*nodes-size(unknown_u_index,1))';    %u_known为已知的节点位移。


%求解插值函数的基本思路：
%   首先利用各节点的自然坐标值构建矩阵A，每行是[1,si,ti,zii,si*ti,ti*zii,zii*si,si*ti*zii]。 
%   解方程组A*[a1;a2;a3;a4...a8]=[u1;u2;u3;u4...u8],ui列阵是各节点虚设的u位移，只是一个为计算方便引入的符号。
%   场变量u位移表示为[a1,a2,a3,a4...a8]*[1,s,t,zi,s*t,t*zi,zi*s,s*t*zi]'，从中再整理符号量u1\u2\u3\u4...u8的系数，他们的系数就是Ni。
%   全局坐标x用[Ni]*[x1;x2;x3;x4;...x8]表示，xi是已知的坐标数值，x\y\z就用s\t\zi表示了。
%   以上只是一个单元的情况。这个例子中取的六面体元，各单元尺寸完全一样，各单元对应的[Ni]是一样的。
%   不规则情况，各单元的节点值只要还是和六面体元一样自然坐标按八个±1取，[Ni]也是一样的。


syms s t zi u1 u2 u3 u4 u5 u6 u7 u8
x_expression=sym(zeros(eles,1));     %x符号列阵用来存储各单元的场变量x以s\t\zi表示的表达式。y符号矩阵对应。
y_expression=sym(zeros(eles,1));
z_expression=sym(zeros(eles,1));


stzi_coordinate=sym([1,1,1;-1,1,1;-1,-1,1;1,-1,1;1,1,-1;-1,1,-1;-1,-1,-1;1,-1,-1]');       %stzi_coordinate符号列阵表示在自然坐标下，每个单元对应的节点的坐标。
u=sym([u1;u2;u3;u4;u5;u6;u7;u8]);     %u以符号表示单元四个节点的u位移，只作为符号在计算中用来求插值函数，没有实际对应的已知值。


stzi_ele_mat=sym(zeros(8,8));         %st_ele_mat是在A*[a1;a2;a3;a4..a8]=[u1;u2;u3;u4...u8]方程组中的A，每行是[1,si,ti,zii,si*ti,ti*zii,zii*si,si*ti*zii]。
for i=1:8
    si=stzi_coordinate(1,i);
    ti=stzi_coordinate(2,i);
    zii=stzi_coordinate(3,i);
    stzi_ele_mat(i,:)=[1,si,ti,zii,si*ti,ti*zii,zii*si,si*ti*zii];
end

 
ai=stzi_ele_mat\u;        %ai为[a1;a2;a3;a4...a8]矩阵。
u_expression=ai'*sym([1,s,t,zi,s*t,t*zi,zi*s,s*t*zi]');       %u_expression是以s\t\zi表示的场变量u函数。
Ni=brick_coeffs_ui(u_expression,u);  %生成Ni，依靠对ui的系数的提取。

x_ele=sym(zeros(8,eles));
y_ele=sym(zeros(8,eles));   %x_ele是各单元节点的总体坐标x的集合，是已知的数值。
z_ele=sym(zeros(8,eles));
for i=1:eles
    x_ele(:,i)=node_coordinate([1,4,7,10,13,16,19,22],i);
    y_ele(:,i)=node_coordinate([2,5,8,11,14,17,20,23],i);
    z_ele(:,i)=node_coordinate([3,6,9,12,15,18,21,24],i);
end


for i=1:eles
    x_expression(i)=Ni'*x_ele(:,i);
    y_expression(i)=Ni'*y_ele(:,i);
    z_expression(i)=Ni'*z_ele(:,i);
end     %x、y的s\t表达式已经得到。


b=sym(zeros(6,24,eles));     %b为各单元B矩阵的集合。
for i=1:eles
    for j=1:8
        b(:,3*j-2:3*j,i)=brick_ele_b_1(Ni(j),x_expression(i),y_expression(i),z_expression(i));
    end
end

 
d=e/((1+nu)*(1-2*nu))*[1-nu,nu,nu,0,0,0;nu,1-nu,nu,0,0,0;nu,nu,1-nu,0,0,0;
     0,0,0,(1-2*nu)/2,0,0;0,0,0,0,(1-2*nu)/2,0;0,0,0,0,0,(1-2*nu)/2];
 
 kk=sym(zeros(24,24,eles));    %kk是各单元刚度矩阵的集合。
 for i=1:1:eles
      kk(:,:,i)=brick_ele_stiff_mat_1(b(:,:,i),d,x_expression(i),y_expression(i),z_expression(i));
 end
kk=double(kk);
 
k=zeros(3*nodes,3*nodes);   %k为总体刚度矩阵。
for i=1:1:eles
    k=brick_assembly_mat(k,kk(:,:,i),node_number(:,i));
end


unknown_var=tetra_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
ua = unknown_var{1,1};
ua=double(ua);
fc = unknown_var{1,2};
fc=double(fc);
U = unknown_var{1,3};
F = unknown_var{1,4};
kcc = unknown_var{1,5};
kca = unknown_var{1,6};
kac = unknown_var{1,7};
kaa = unknown_var{1,8};
 
u=zeros(24,eles);    %u为各单元对应节点的位移列阵集合。
for i=1:1:eles
      u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3);U(node_number(3,i)*3-2:node_number(3,i)*3);
        U(node_number(4,i)*3-2:node_number(4,i)*3);U(node_number(5,i)*3-2:node_number(5,i)*3);
        U(node_number(6,i)*3-2:node_number(6,i)*3);U(node_number(7,i)*3-2:node_number(7,i)*3);U(node_number(8,i)*3-2:node_number(8,i)*3)];
end
 
 stress=sym(zeros(6,eles));   %sigma为单元应力矩阵，此时算出的是关于坐标的变量，不再是常量，strain为应变。
 strain=sym(zeros(6,eles));
 for i=1:1:eles
     stress(:,i)=d*b(:,:,i)*u(:,i);
     strain(:,i)=b(:,:,i)*u(:,i);
 end
 num_stress=vpa(stress,6);
 num_strain=vpa(strain,6);
 stress_zero_point=double(subs(stress,{s,t,zi},{0,0,0}));
 strain_zero_point=double(subs(strain,{s,t,zi},{0,0,0}));

 principal_stress=zeros(3,eles);    %求解主应力与主方向。
 principal_strain=zeros(3,eles);
 principal_direction=zeros(3,3,eles);
  for i=1:1:eles
      m_s=tetra_main_stress(stress_zero_point(:,i));
      m_strain=tetra_main_stress(strain_zero_point(:,i));
      principal_stress(:,i)=m_s{1,1};
      principal_direction(:,:,i)=m_s{1,2};
      principal_strain(:,i)=m_strain{1,1};
  end
 


%以下写文件。
f=fopen('brick_eg_output.xls','w');

fprintf(f, 'element stiffness matrices:\r\n\r\n');
for i=1:eles
    for j=1:24
        for m=1:24
            fprintf(f,'%8.4f\t',kk(j,m,i));
        end
        fprintf(f, '\r\n');
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\nglobal stiffness matrix:\r\n\r\n');
for j=1:3*nodes
    for m=1:3*nodes
        fprintf(f,'%8.4f\t',k(j,m));
    end
    fprintf(f, '\r\n');
end

fprintf(f,'\r\n\r\nK_cc= \r\n\r\n');
for i=1:size(kcc,1)
    for j=1:size(kcc,2)
        fprintf(f,'%8.4f\t',kcc(i,j));
    end
    fprintf(f,'\r\n');
end
fprintf(f,'\r\n\r\nK_ca= \r\n\r\n');
for i=1:size(kca,1)
    for j=1:size(kca,2)
        fprintf(f,'%8.4f\t',kca(i,j));
    end
    fprintf(f,'\r\n');
end
fprintf(f,'\r\n\r\nK_ac= \r\n\r\n');
for i=1:size(kac,1)
    for j=1:size(kac,2)
        fprintf(f,'%8.4f\t',kac(i,j));
    end
    fprintf(f,'\r\n');
end
fprintf(f,'\r\n\r\nK_aa= \r\n\r\n');
for i=1:size(kaa,1)
    for j=1:size(kaa,2)
        fprintf(f,'%8.4f\t',kaa(i,j));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nU_a= \r\n\r\n');
for i=1:size(ua,1)
    fprintf(f,'%8.4f\t',ua(i));
end
fprintf(f,'\r\n\r\nF_c= \r\n\r\n');
for i=1:size(fc,1)
    fprintf(f,'%8.4f\t',fc(i));
end
fprintf(f,'\r\n\r\nU= \r\n\r\n');
for i=1:size(U,1)
    fprintf(f,'%8.4f\t',U(i));
end
fprintf(f,'\r\n\r\nF= \r\n\r\n');
for i=1:size(F,1)
    fprintf(f,'%8.4f\t',F(i));
end

fprintf(f,'\r\n\r\nstress= \r\n\r\n');
for i=1:size(num_stress,1)
    for j=1:size(num_stress,2)
        fprintf(f,'%s\t\t',char(num_stress(i,j)));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nstrain= \r\n\r\n');
for i=1:size(num_strain,1)
    for j=1:size(num_strain,2)
        fprintf(f,'%s\t\t',char(num_strain(i,j)));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nstress_zero_point= \r\n\r\n');
for i=1:size(stress_zero_point,1)
    for j=1:size(stress_zero_point,2)
        fprintf(f,'%8.4f\t\t',stress_zero_point(i,j));
    end
    fprintf(f,'\r\n');
end


fprintf(f,'\r\n\r\nstrain_zero_point= \r\n\r\n');
for i=1:size(strain_zero_point,1)
    for j=1:size(strain_zero_point,2)
        fprintf(f,'%8.4f\t\t',strain_zero_point(i,j));
    end
    fprintf(f,'\r\n');
end


fprintf(f,'\r\n\r\nprincipal_stress= \r\n\r\n');
for i=1:size(principal_stress,1)
    for j=1:size(principal_stress,2)
        fprintf(f,'%8.4f\t\t',principal_stress(i,j));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nprincipal_strain= \r\n\r\n');
for i=1:size(principal_strain,1)
    for j=1:size(principal_strain,2)
        fprintf(f,'%8.4f\t\t',principal_strain(i,j));
    end
    fprintf(f,'\r\n');
end


fprintf(f,'\r\n\r\nprincipal_direction= \r\n\r\n');
for p=1:eles
    for i=1:size(principal_direction,1)
        for j=1:size(principal_direction,2)
            fprintf(f,'%8.4f\t\t',principal_direction(i,j,p));
        end
        fprintf(f,'\r\n');
    end
    fprintf(f,'\r\n');
end
fclose(f);




% e=210e9;
% nu=0.3;
% eles=2;
% nodes=12;
% node_number=[11,12,5,6,9,10,1,2;
%     8,7,12,11,4,3,10,9]';  %每个单元对应节点的编号，按一到八卦限排列。
% node_coordinate=[0,0,0;
%     0.025,0,0;
%     0,0.5,0;
%     0.025,0.5,0;
%     0,0,0.25;
%     0.025,0,0.25;
%     0,0.5,0.25;
%     0.025,0.5,0.25;
%     0.025,0.25,0;
%     0,0.25,0;
%     0.025,0.25,0.25;
%     0,0.25,0.25]';  %各节点的坐标列阵。
% node_coordinate=brick_ele_coordinate(node_coordinate,node_number);  %生成每个单元对应节点的坐标列阵集合。
% syms s t zi
% s_nodes=[1;-1;-1;1;1;-1;-1;1];   %s_nodes针对局部坐标系下一个单元各节点的横坐标，每个单元都一样,是常值。t_nodes对应纵坐标，zi对应竖坐标。
% t_nodes=[1;1;-1;-1;1;1;-1;-1];
% zi_nodes=[1;1;1;1;-1;-1;-1;-1];
% n=sym(zeros(8,eles));   %n为各单元插值函数的集合，sym类型内容，8行eles列。
% s_eles=sym(zeros(eles,1)); %xi_eles针对各个单元中的无量纲横坐标x/a，每个单元因为a可能不一样，各元素也不一定一样。sym类型内容。eta_eles对应纵坐标。
% t_eles=sym(zeros(eles,1));
% zi_eles=sym(zeros(eles,1));
% for i=1:1:eles
%     s_eles(i)=s;
%     t_eles(i)=t;
%     zi_eles(i)=zi;
% end
% for i=1:1:8
%     for j=1:1:eles
%     n(i,j)=0.125*(1+s_nodes(i)*s_eles(j))*(1+t_nodes(i)*t_eles(j))*(1+zi_nodes(i)*zi_eles(j));
%     end
% end
% jacob_eles=sym(zeros(8,3*eles));    %生成各单元对应的插值函数组的雅可比矩阵集合。
% for i=1:1:eles
%     jacob_eles(:,3*i-2:3*i)=jacobian(n(:,i));
% end
% 
% x=sym(zeros(eles,1));  %x是各单元插值后以自然坐标系变量描述的总体坐标系x变量。
% y=sym(zeros(eles,1));
% z=sym(zeros(eles,1));
% for i=1:1:eles
%     xyz=brick_xyz_eles(n(:,i),node_coordinate(:,i));
%     x(i)=xyz(1);
%     y(i)=xyz(2);
%     z(i)=xyz(3);
% end
% 
% b=sym(zeros(6*eles,24));     %b为各单元B矩阵的集合。
% for i=1:1:eles
%     for j=1:1:8
%     b(6*i-5:6*i,3*j-2:3*j)=brick_ele_b(n(j,i),x(i),y(i),z(i));
%     end
% end     %生成B矩阵集合。
% 
% d=e/((1+nu)*(1-2*nu))*[1-nu,nu,nu,0,0,0;nu,1-nu,nu,0,0,0;nu,nu,1-nu,0,0,0;
%     0,0,0,(1-2*nu)/2,0,0;0,0,0,0,(1-2*nu)/2,0;0,0,0,0,0,(1-2*nu)/2];
% 
% kk=sym(zeros(24*eles,24));    %kk是各单元刚度矩阵的集合。
% for i=1:1:eles
%      kk(24*i-23:24*i,:)=brick_ele_stiff_mat(b(6*i-5:6*i,:),d,x(i),y(i),z(i));
% end
% kk=double(kk)
% 
%  k=zeros(3*nodes,3*nodes);   %k为总体刚度矩阵。
%  for i=1:1:eles
%      k=brick_assembly_mat(k,kk(24*i-23:24*i,:),node_number(:,i));
%  end
%  unknown_u_index=[7,8,9,10,11,12,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]';   %unknown_u_index是未知位移序号的索引。
%   known_f=[0;4.6875e3;0;0;4.6875e3;0;0;4.6875e3;0;0;4.6875e3;0;zeros(12,1)];
%   known_u=zeros(12,1);
%   unknown_var=tetra_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
%   unknown_u=unknown_var{1,1}
%   unknown_f=unknown_var{1,2}
%  U=unknown_var{1,3}  %U为总位移列阵。
%  F=unknown_var{1,4}   %F为总受力列阵。
 
%  u=zeros(24,eles);    %u为各单元对应节点的位移列阵集合。
%  for i=1:1:eles
%        u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3);U(node_number(3,i)*3-2:node_number(3,i)*3);
%          U(node_number(4,i)*3-2:node_number(4,i)*3);U(node_number(5,i)*3-2:node_number(5,i)*3);
%          U(node_number(6,i)*3-2:node_number(6,i)*3);U(node_number(7,i)*3-2:node_number(7,i)*3);U(node_number(8,i)*3-2:node_number(8,i)*3)];
%  end
%  
%  stress=sym(zeros(6,eles));   %sigma为单元应力矩阵，此时算出的是关于坐标的变量，不再是常量，strain为应变。
%  strain=sym(zeros(6,eles));
%  for i=1:1:eles
%      stress(:,i)=d*b(6*i-5:6*i,:)*u(:,i);
%      strain(:,i)=b(6*i-5:6*i,:)*u(:,i);
%  end
%  num_stress=vpa(stress,6)
%  num_strain=vpa(strain,6)
%  stress_zero_point=double(subs(stress,{s,t,zi},{0,0,0}))
%  strain_zero_point=double(subs(strain,{s,t,zi},{0,0,0}))
% 
%  principal_stress=zeros(3,eles);    %求解主应力与主方向。
%  principal_strain=zeros(3,eles);
%  principal_direction=zeros(3*eles,3);
%   for i=1:1:eles
%       m_s=tetra_main_stress(stress_zero_point(:,i));
%       m_strain=tetra_main_stress(strain_zero_point(:,i));
%       principal_stress(:,i)=m_s{1,1};
%       principal_direction(3*i-2:3*i,:)=m_s{1,2};
%       principal_strain(:,i)=m_strain{1,1};
%   end
%   principal_stress
%   principal_strain
%   principal_direction
%  
%  
 
 
