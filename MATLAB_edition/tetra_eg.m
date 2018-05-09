%������Ԫ������

source=xlsread('tetra_eg_input.xlsx');
source(all(isnan(source),2),:)=[];
eles=source(1,1);
nodes=source(1,2);
e=source(2,1);
nu=source(2,2);
node_number=source(3:2+eles,1:4)';   %ÿ����Ԫ��Ӧ�ڵ�ı�š�
node_coordinate=source(3+eles:2+eles+nodes,1:3)';
node_coordinate=tetra_ele_coordinate(node_coordinate,node_number);   %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
unknown_u_index=source(3+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';       %unknown_u_index��U��������δ֪λ��λ�õ�������
known_f=source(4+eles+nodes,1:size(unknown_u_index,1))';     %f_knownΪ��֪�Ľڵ�����������Ч�ڵ��غɡ�
known_u=source(5+eles+nodes,1:3*nodes-size(unknown_u_index,1))';    %u_knownΪ��֪�Ľڵ�λ�ơ�

alpha=zeros(4,eles);    %�����ֵ��������Ҫϵ��������
beta=zeros(4,eles);
gamma=zeros(4,eles);
delta=zeros(4,eles);
for i=1:1:eles
    alpha(:,i)=tetra_alpha(alpha(:,i),node_coordinate(:,i));
    beta(:,i)=tetra_beta(beta(:,i),node_coordinate(:,i));
    gamma(:,i)=tetra_gamma(gamma(:,i),node_coordinate(:,i));
    delta(:,i)=tetra_delta(delta(:,i),node_coordinate(:,i));
end

v=zeros(1,eles);    %�����Ԫ�������
for i=1:1:eles
    v(i)=tetra_volume(node_coordinate(:,i));
end

b=zeros(6,12,eles);  %bΪ����ԪB����ļ��ϡ�
for i=1:1:eles
    for j=1:1:4
    b(:,3*j-2:3*j,i)=tetra_ele_b(v(i),beta(j,i),gamma(j,i),delta(j,i));
    end
end     %���ɸ���ԪB����ļ��ϡ�

d=e/((1+nu)*(1-2*nu))*[1-nu,nu,nu,0,0,0;nu,1-nu,nu,0,0,0;nu,nu,1-nu,0,0,0;
    0,0,0,(1-2*nu)/2,0,0;0,0,0,0,(1-2*nu)/2,0;0,0,0,0,0,(1-2*nu)/2];
kk=zeros(12,12,eles);     %kkΪ����Ԫ�նȾ���ļ��ϡ�
for i=1:1:eles
    kk(:,:,i)=tetra_ele_stiff_mat(b(:,:,i),d,v(i));
end     %��ɵ�Ԫ�նȾ��󼯺ϡ�


k=zeros(3*nodes,3*nodes);   %kΪ����նȾ���
for i=1:1:eles
    k=tetra_assembly_mat(k,kk(:,:,i),node_number(:,i));
end     %�������նȾ���


unknown_var=tetra_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
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

u=zeros(12,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
for i=1:1:eles
     u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3);U(node_number(3,i)*3-2:node_number(3,i)*3);
         U(node_number(4,i)*3-2:node_number(4,i)*3)];
end
   
stress=zeros(6,eles);    %sigmaΪÿ��Ԫ��Ӧ�����С�
strain=zeros(6,eles);
for i=1:1:eles
   stress(:,i)=d*b(:,:,i)*u(:,i);
   strain(:,i)=b(:,:,i)*u(:,i);
end

principal_stress=zeros(3,eles);    %�����Ӧ����������
principal_strain=zeros(3,eles);
principal_direction=zeros(3,3,eles);
for i=1:1:eles
    m_s=tetra_main_stress(stress(:,i));
    m_strain=tetra_main_stress(strain(:,i));
    principal_stress(:,i)=m_s{1,1};
    principal_direction(:,:,i)=m_s{1,2};
    principal_strain(:,i)=m_strain{1,1};
end

%����д�ļ���
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
%     1,6,4,7]';  %ÿ����Ԫ��Ӧ�ڵ�ı�ţ������һ���ڵ㿴ȥ��ǰ�����ڵ㰴��ʱ�����У���֤�������ʽֵΪ����
% node_coordinate=[0,0,0;
%     0.025,0,0;
%     0,0.5,0;
%     0.025,0.5,0;
%     0,0,0.25;
%     0.025,0,0.25;
%     0,0.5,0.25;
%     0.025,0.5,0.25]';   %ÿ���ڵ��Ӧ�ı�š�
% node_coordinate=tetra_ele_coordinate(node_coordinate,node_number);  %����ÿ����Ԫ��Ӧ�ڵ���������󼯺ϡ�
% alpha=zeros(4,eles);    %�����ֵ��������Ҫϵ��������
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
% v=zeros(1,eles);    %�����Ԫ�������
% for i=1:1:eles
%     v(i)=tetra_volume(node_coordinate(:,i));
% end
% 
% b=zeros(6*eles,12);  %bΪ����ԪB����ļ��ϡ�
% for i=1:1:eles
%     for j=1:1:4
%     b(6*i-5:6*i,3*j-2:3*j)=tetra_ele_b(v(i),beta(j,i),gamma(j,i),delta(j,i));
%     end
% end     %���ɸ���ԪB����ļ��ϡ�
% 
% d=e/((1+nu)*(1-2*nu))*[1-nu,nu,nu,0,0,0;nu,1-nu,nu,0,0,0;nu,nu,1-nu,0,0,0;
%     0,0,0,(1-2*nu)/2,0,0;0,0,0,0,(1-2*nu)/2,0;0,0,0,0,0,(1-2*nu)/2];
% kk=zeros(12*eles,12);     %kkΪ����Ԫ�նȾ���ļ��ϡ�
% for i=1:1:eles
%     kk(12*i-11:12*i,:)=tetra_ele_stiff_mat(b(6*i-5:6*i,:),d,v(i));
% end     %��ɵ�Ԫ�նȾ��󼯺ϡ�
% 
% k=zeros(3*nodes,3*nodes);   %kΪ����նȾ���
% for i=1:1:eles
%     k=tetra_assembly_mat(k,kk(12*i-11:12*i,:),node_number(:,i));
% end     %�������նȾ���
%  unknown_u_index=[7,8,9,10,11,12,19,20,21,22,23,24]';   %unknown_u_index��δ֪λ����ŵ�������
%  known_f=[0;3.125e3;0;0;6.25e3;0;0;6.25e3;0;0;3.125e3;0];
%  known_u=zeros(12,1);
%  unknown_var=tetra_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
%  unknown_u=unknown_var{1,1}
%  unknown_f=unknown_var{1,2}
% U=unknown_var{1,3}  %UΪ��λ������
% F=unknown_var{1,4}   %FΪ����������
% u=zeros(12,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
% 
% for i=1:1:eles
%      u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3);U(node_number(3,i)*3-2:node_number(3,i)*3);
%          U(node_number(4,i)*3-2:node_number(4,i)*3)];
% end
%   
% sigma=zeros(6,eles);    %sigmaΪÿ��Ԫ��Ӧ�����С�
% strain=zeros(6,eles);
%  for i=1:1:eles
%     sigma(:,i)=d*b(6*i-5:6*i,:)*u(:,i);
%     strain(:,i)=b(6*i-5:6*i,:)*u(:,i);
%  end
%  sigma
%  strain
%   principal_stress=zeros(3,eles);    %�����Ӧ����������
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
