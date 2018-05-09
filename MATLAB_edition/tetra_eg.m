%四面体元算例。

source=xlsread('tetra_eg_input.xlsx');
source(all(isnan(source),2),:)=[];
eles=source(1,1);
nodes=source(1,2);
e=source(2,1);
nu=source(2,2);
node_number=source(3:2+eles,1:4)';   %每个单元对应节点的编号。
node_coordinate=source(3+eles:2+eles+nodes,1:3)';
node_coordinate=tetra_ele_coordinate(node_coordinate,node_number);   %每个单元对应节点的坐标。
unknown_u_index=source(3+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';       %unknown_u_index是U列向量中未知位移位置的索引。
known_f=source(4+eles+nodes,1:size(unknown_u_index,1))';     %f_known为已知的节点力，包含等效节点载荷。
known_u=source(5+eles+nodes,1:3*nodes-size(unknown_u_index,1))';    %u_known为已知的节点位移。

alpha=zeros(4,eles);    %定义插值函数中重要系数参量。
beta=zeros(4,eles);
gamma=zeros(4,eles);
delta=zeros(4,eles);
for i=1:1:eles
    alpha(:,i)=tetra_alpha(alpha(:,i),node_coordinate(:,i));
    beta(:,i)=tetra_beta(beta(:,i),node_coordinate(:,i));
    gamma(:,i)=tetra_gamma(gamma(:,i),node_coordinate(:,i));
    delta(:,i)=tetra_delta(delta(:,i),node_coordinate(:,i));
end

v=zeros(1,eles);    %求各单元的体积。
for i=1:1:eles
    v(i)=tetra_volume(node_coordinate(:,i));
end

b=zeros(6,12,eles);  %b为各单元B矩阵的集合。
for i=1:1:eles
    for j=1:1:4
    b(:,3*j-2:3*j,i)=tetra_ele_b(v(i),beta(j,i),gamma(j,i),delta(j,i));
    end
end     %生成各单元B矩阵的集合。

d=e/((1+nu)*(1-2*nu))*[1-nu,nu,nu,0,0,0;nu,1-nu,nu,0,0,0;nu,nu,1-nu,0,0,0;
    0,0,0,(1-2*nu)/2,0,0;0,0,0,0,(1-2*nu)/2,0;0,0,0,0,0,(1-2*nu)/2];
kk=zeros(12,12,eles);     %kk为各单元刚度矩阵的集合。
for i=1:1:eles
    kk(:,:,i)=tetra_ele_stiff_mat(b(:,:,i),d,v(i));
end     %完成单元刚度矩阵集合。


k=zeros(3*nodes,3*nodes);   %k为总体刚度矩阵。
for i=1:1:eles
    k=tetra_assembly_mat(k,kk(:,:,i),node_number(:,i));
end     %完成总体刚度矩阵。


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

u=zeros(12,eles);    %u为各单元对应节点的位移列阵集合。
for i=1:1:eles
     u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3);U(node_number(3,i)*3-2:node_number(3,i)*3);
         U(node_number(4,i)*3-2:node_number(4,i)*3)];
end
   
stress=zeros(6,eles);    %sigma为每单元的应力阵列。
strain=zeros(6,eles);
for i=1:1:eles
   stress(:,i)=d*b(:,:,i)*u(:,i);
   strain(:,i)=b(:,:,i)*u(:,i);
end

principal_stress=zeros(3,eles);    %求解主应力与主方向。
principal_strain=zeros(3,eles);
principal_direction=zeros(3,3,eles);
for i=1:1:eles
    m_s=tetra_main_stress(stress(:,i));
    m_strain=tetra_main_stress(strain(:,i));
    principal_stress(:,i)=m_s{1,1};
    principal_direction(:,:,i)=m_s{1,2};
    principal_strain(:,i)=m_strain{1,1};
end

%以下写文件。
f=fopen('tetra_eg_output.txt','w');

fprintf(f, 'element stiffness matrices:\r\n\r\n');
for i=1:eles
    for j=1:12
        for m=1:12
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
for i=1:size(stress,1)
    for j=1:size(stress,2)
        fprintf(f,'%8.4f\t\t',stress(i,j));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nstrain= \r\n\r\n');
for i=1:size(strain,1)
    for j=1:size(strain,2)
        fprintf(f,'%8.4f\t\t',strain(i,j));
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
    for i=1:size(principal_direction(:,:,p),1)
        for j=1:size(principal_direction(:,:,p),2)
            fprintf(f,'%8.4f\t',principal_direction(i,j,p));
        end
        fprintf(f,'\r\n');
    end
    fprintf(f,'\r\n');
end
fclose(f);







% e=210e9;
% nu=0.3;
% eles=5;
% nodes=8;
% node_number=[1,2,4,6;
%     1,4,3,7;
%     6,5,7,1;
%     6,7,8,4;
%     1,6,4,7]';  %每个单元对应节点的编号，从最后一个节点看去，前三个节点按逆时针排列，保证体积行列式值为正。
% node_coordinate=[0,0,0;
%     0.025,0,0;
%     0,0.5,0;
%     0.025,0.5,0;
%     0,0,0.25;
%     0.025,0,0.25;
%     0,0.5,0.25;
%     0.025,0.5,0.25]';   %每个节点对应的编号。
% node_coordinate=tetra_ele_coordinate(node_coordinate,node_number);  %生成每个单元对应节点的坐标列阵集合。
% alpha=zeros(4,eles);    %定义插值函数中重要系数参量。
% beta=zeros(4,eles);
% gamma=zeros(4,eles);
% delta=zeros(4,eles);
% for i=1:1:eles
%     alpha(:,i)=tetra_alpha(alpha(:,i),node_coordinate(:,i));
%     beta(:,i)=tetra_beta(beta(:,i),node_coordinate(:,i));
%     gamma(:,i)=tetra_gamma(gamma(:,i),node_coordinate(:,i));
%     delta(:,i)=tetra_delta(delta(:,i),node_coordinate(:,i));
% end
% 
% v=zeros(1,eles);    %求各单元的体积。
% for i=1:1:eles
%     v(i)=tetra_volume(node_coordinate(:,i));
% end
% 
% b=zeros(6*eles,12);  %b为各单元B矩阵的集合。
% for i=1:1:eles
%     for j=1:1:4
%     b(6*i-5:6*i,3*j-2:3*j)=tetra_ele_b(v(i),beta(j,i),gamma(j,i),delta(j,i));
%     end
% end     %生成各单元B矩阵的集合。
% 
% d=e/((1+nu)*(1-2*nu))*[1-nu,nu,nu,0,0,0;nu,1-nu,nu,0,0,0;nu,nu,1-nu,0,0,0;
%     0,0,0,(1-2*nu)/2,0,0;0,0,0,0,(1-2*nu)/2,0;0,0,0,0,0,(1-2*nu)/2];
% kk=zeros(12*eles,12);     %kk为各单元刚度矩阵的集合。
% for i=1:1:eles
%     kk(12*i-11:12*i,:)=tetra_ele_stiff_mat(b(6*i-5:6*i,:),d,v(i));
% end     %完成单元刚度矩阵集合。
% 
% k=zeros(3*nodes,3*nodes);   %k为总体刚度矩阵。
% for i=1:1:eles
%     k=tetra_assembly_mat(k,kk(12*i-11:12*i,:),node_number(:,i));
% end     %完成总体刚度矩阵。
%  unknown_u_index=[7,8,9,10,11,12,19,20,21,22,23,24]';   %unknown_u_index是未知位移序号的索引。
%  known_f=[0;3.125e3;0;0;6.25e3;0;0;6.25e3;0;0;3.125e3;0];
%  known_u=zeros(12,1);
%  unknown_var=tetra_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
%  unknown_u=unknown_var{1,1}
%  unknown_f=unknown_var{1,2}
% U=unknown_var{1,3}  %U为总位移列阵。
% F=unknown_var{1,4}   %F为总受力列阵。
% u=zeros(12,eles);    %u为各单元对应节点的位移列阵集合。
% 
% for i=1:1:eles
%      u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3);U(node_number(3,i)*3-2:node_number(3,i)*3);
%          U(node_number(4,i)*3-2:node_number(4,i)*3)];
% end
%   
% sigma=zeros(6,eles);    %sigma为每单元的应力阵列。
% strain=zeros(6,eles);
%  for i=1:1:eles
%     sigma(:,i)=d*b(6*i-5:6*i,:)*u(:,i);
%     strain(:,i)=b(6*i-5:6*i,:)*u(:,i);
%  end
%  sigma
%  strain
%   principal_stress=zeros(3,eles);    %求解主应力与主方向。
%   principal_strain=zeros(3,eles);
%   principal_direction=zeros(3*eles,3);
%   for i=1:1:eles
%       m_s=tetra_main_stress(sigma(:,i));
%       m_strain=tetra_main_stress(strain(:,i));
%       principal_stress(:,i)=m_s{1,1};
%       principal_direction(3*i-2:3*i,:)=m_s{1,2};
%       principal_strain(:,i)=m_strain{1,1};
%   end
%   principal_stress
%   principal_strain
%   principal_direction
