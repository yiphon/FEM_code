%������Ԫ������

s=xlsread('cst_eg_input.xlsx');
s(all(isnan(s),2),:)=[];
eles=s(1,1);
nodes=s(1,2);
e=s(2,1);
nu=s(2,2);
t=s(2,3);
node_number=s(3:2+eles,1:3)';   %ÿ����Ԫ��Ӧ�ڵ�ı�š�
node_coordinate=s(3+eles:2+eles+nodes,1:2)';
node_coordinate=tri_ele_coordinate(node_coordinate,node_number);   %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
unknown_u_index=s(3+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';       %unknown_u_index��U��������δ֪λ��λ�õ�������
f_known=s(4+eles+nodes,1:size(unknown_u_index,1))';     %f_knownΪ��֪�Ľڵ�����������Ч�ڵ��غɡ�
u_known=s(5+eles+nodes,1:2*nodes-size(unknown_u_index,1))';    %u_knownΪ��֪�Ľڵ�λ�ơ�

a=zeros(eles,1);
for i=1:1:eles
    a(i)=tri_ele_square(node_coordinate(:,i));
end     %���ɸ���Ԫ�������

b=zeros(3,6,eles);  %bΪ����ԪB����ļ��ϡ�
for i=1:1:eles
    b(:,:,i)=tri_ele_b(a(i),node_coordinate(:,i));
end     %���ɸ���ԪB����ļ��ϡ�

d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
kk=zeros(6,6,eles);     %kkΪ����Ԫ�նȾ���ļ��ϡ�
for i=1:1:eles
    kk(:,:,i)=tri_ele_stiff_mat(a(i),t,b(:,:,i),d);
end     %��ɵ�Ԫ�նȾ��󼯺ϡ�

k=zeros(2*nodes,2*nodes);   %kΪ����նȾ���
for i=1:1:eles
    k=tri_assembly_mat(k,kk(:,:,i),node_number(:,i));
end     %�������նȾ���

unknown_var=tri_solve_unknown(k,unknown_u_index,f_known,u_known,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
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

u=zeros(6,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
for i=1:1:eles
     u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2)];
end

sigma=zeros(3,eles);    %sigmaΪÿ��Ԫ��Ӧ�����С�
 for i=1:1:eles
    sigma(:,i)=d*b(:,:,i)*u(:,i);
 end
 
 main_stress=zeros(3,eles);    %�����Ӧ����������
 for i=1:1:eles
     main_stress(:,i)=tri_main_sigma(sigma(:,i));
 end
principal_stress_=main_stress(1:2,:);
principal_direction=main_stress(3,:);


%����д�ļ���
f=fopen('cst_eg_output.xls','w');

fprintf(f, 'element stiffness matrices:\r\n\r\n');
for i=1:eles
    for j=1:6
        for m=1:6
            fprintf(f,'%8.4f\t',kk(j,m,i));
        end
        fprintf(f, '\r\n');
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\nglobal stiffness matrix:\r\n\r\n');
for j=1:2*nodes
    for m=1:2*nodes
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

fprintf(f,'\r\n\r\nsigma= \r\n\r\n');
for i=1:size(sigma,1)
    for j=1:size(sigma,2)
        fprintf(f,'%8.4f\t\t',sigma(i,j));
    end
    fprintf(f,'\r\n');
end


fprintf(f,'\r\n\r\nprincipal_stress= \r\n\r\n');
for i=1:size(principal_stress,1)-1
    for j=1:size(principal_stress,2)
        fprintf(f,'%8.4f\t\t',principal_stress(i,j));
    end
    fprintf(f,'\r\n');
end

fprintf(f,'\r\n\r\nprincipal_direction= \r\n\r\n');
for i=1:size(principal_direction,1)
    for j=1:size(principal_direction,2)
        fprintf(f,'%8.4f\t\t',principal_direction(i,j));
    end
    fprintf(f,'\r\n');
end

fclose(f);







% e=210e9;
% nu=0.3;
% t=0.025;
% eles=12;
% nodes=11;
% node_number=[1,3,2;
%     1,4,3;
%     2,3,5;
%     3,4,5;
%     4,6,5;
%     4,7,6;
%     5,6,8;
%     6,7,8;
%     5,9,10;
%     5,8,9;
%     9,11,10;
%     9,8,11]';  %ÿ����Ԫ��Ӧ�ڵ�ı�ţ�����ʱ�����С�
% a1=[0,0.5];
% a2=[0.25,0.5];
% a3=[0.125,0.375];
% a4=[0,0.25];
% a5=[0.25,0.25];
% a6=[0.125,0.125];
% a7=[0,0];
% a8=[0.25,0];
% a9=[0.375,0.125];
% a10=[0.5,0.25];
% a11=[0.5,0];
% node_coordinate=[a1,a3,a2;
%     a1,a4,a3;
%     a2,a3,a5;
%     a3,a4,a5;
%     a4,a6,a5;
%     a4,a7,a6;
%     a5,a6,a8;
%     a6,a7,a8;
%     a5,a9,a10;
%     a5,a8,a9;
%     a9,a11,a10;
%     a9,a8,a11;]';  %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
% a=zeros(eles,1);
% for i=1:1:eles
%     a(i)=tri_ele_square(node_coordinate(:,i));
% end     %���ɸ���Ԫ�������
% b=zeros(3*eles,6);  %bΪ����ԪB����ļ��ϡ�
% for i=1:1:eles
%     b(3*i-2:3*i,:)=tri_ele_b(a(i),node_coordinate(:,i));
% end     %���ɸ���ԪB����ļ��ϡ�
% d=e/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
% kk=zeros(6*eles,6);     %kkΪ����Ԫ�նȾ���ļ��ϡ�
% for i=1:1:eles
%     kk(6*i-5:6*i,:)=tri_ele_stiff_mat(a(i),t,b(3*i-2:3*i,:),d);
% end     %��ɵ�Ԫ�նȾ��󼯺ϡ�
% k=zeros(2*nodes,2*nodes);   %kΪ����նȾ���
% for i=1:1:eles
%     k=tri_assembly_mat(k,kk(6*i-5:6*i,:),node_number(:,i));
% end     %�������նȾ���
%  unknown_u_index=[3,4,5,6,9,10,11,12,15,16,17,18,19,20,21,22]';   %unknown_u_index��δ֪λ����ŵ�������
%  known_f=[zeros(5,1);-12.5e3;zeros(10,1)];
%  known_u=zeros(6,1);
%  unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
%  unknown_u=unknown_var{1,1}
%  unknown_f=unknown_var{1,2}
% U=unknown_var{1,3}  %UΪ��λ������
% F=unknown_var{1,4}   %FΪ����������
% u=zeros(6,eles);    %uΪ����Ԫ��Ӧ�ڵ��λ�����󼯺ϡ�
% for i=1:1:eles
%      u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2);U(node_number(3,i)*2-1:node_number(3,i)*2)];
%   end
% sigma=zeros(3,eles);    %sigmaΪÿ��Ԫ��Ӧ�����С�
%  for i=1:1:eles
%     sigma(:,i)=d*b(3*i-2:3*i,:)*u(:,i);
%  end
%  sigma
%  principal_stress=zeros(3,eles);    %�����Ӧ����������
%  for i=1:1:eles
%      principal_stress(:,i)=tri_main_sigma(sigma(:,i));
%  end
%   principal_direction=principal_stress(3,:)
% principal_stress_=principal_stress(1:2,:)
% 








