%����ռ�Ԫ������

s=xlsread('frame3d_eg_input.xlsx');
s(all(isnan(s),2),:)=[];
eles=s(1,1);
nodes=s(1,2);
e=s(2,1:eles);
g=s(3,1:eles);
a=s(4,1:eles);
iy=s(5,1:eles);
iz=s(6,1:eles);
j=s(7,1:eles);
node_number=s(8:7+eles,1:2)';   %ÿ����Ԫ��Ӧ�ڵ�ı�š�
node_coordinate=s(8+eles:7+eles+nodes,1:3)';
node_coordinate=truss3d_ele_coordinate(node_coordinate,node_number);   %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
unknown_u_index=s(8+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';
f_known=s(9+eles+nodes,1:size(unknown_u_index,1))';
u_known=s(10+eles+nodes,1:6*nodes-size(unknown_u_index,1))';

ele_info=zeros(10,eles);    %��Ԫ��Ϣ����ÿ��������e\a\iy\iz\j\l\cXx\cYx\cZx\g��
ele_info(1,:)=e;
ele_info(2,:)=a;
ele_info(3,:)=iy;
ele_info(4,:)=iz;
ele_info(5,:)=j;
ele_info(10,:)=g;
ele_info(6:9,:)=frame3d_ele_info(node_coordinate);  %��ȫ��Ԫ��Ϣ�������Ϣ��

r=zeros(12,12,eles);   %rΪ��Ԫת�ƾ���ļ��ϣ�ÿ��Ԫʮ����12�С�
for i=1:1:eles
    r(:,:,i)=frame3d_r_mat(ele_info(7:9,i));
end     %��ɵ�Ԫת�ƾ���ļ��ϡ�

k_ele=zeros(12,12,eles);  %k_ele�ǵ�Ԫ�ھֲ�����ϵ�µĸնȾ���ÿ��Ԫʮ����12�С�
for i=1:1:eles
    k_ele(:,:,i)=frame3d_k_ele_mat(ele_info(:,i));
end     %��ɾֲ�����ϵ�µĵ�Ԫ�նȾ��󼯺ϡ�

kk=zeros(12,12,eles); %kk�����������µ�Ԫ�նȾ���ļ��ϣ�ÿ��Ԫ12��12�С�
for i=1:1:eles
    kk(:,:,i)=frame3d_ele_stiff_mat(r(:,:,i),k_ele(:,:,i));
end     %�����������ϵ�µ�Ԫ�նȾ��󼯺ϡ�


k=zeros(nodes*6,nodes*6);   %k������նȾ���
for i=1:1:eles
    k=frame3d_assembly_mat(k,kk(:,:,i),node_number(:,i));
end     %w�������նȾ���

unknown_var=frame3d_solve_unknown(k,unknown_u_index,f_known,u_known,nodes);     %���÷ֿ�������δ֪�ı���������λ�ơ�
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

u=zeros(12,eles);   %uΪÿ��Ԫ���ڵ��λ�����С�
for i=1:1:eles
   u(:,i)=[U(node_number(1,i)*6-5:node_number(1,i)*6);U(node_number(2,i)*6-5:node_number(2,i)*6)];
end

f_nodal=zeros(12,eles);    %fΪÿ��Ԫ���ڵ���������У������¶�ÿ����Ԫ�����Ĳ���󣬻�Ӧ��ȥ��Ч�ڵ��غɡ�
for i=1:1:eles
    f_nodal(:,i)=k_ele(:,:,i)*r(:,:,i)*u(:,i);
end


%����д�ļ���
f=fopen('frame3d_eg_output.xls','w');

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
for j=1:6*nodes
    for m=1:6*nodes
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

fprintf(f,'\r\n\r\nf_nodal= \r\n\r\n');
for i=1:size(f_nodal,1)
    for j=1:size(f_nodal,2)
        fprintf(f,'%8.4f\t\t',f_nodal(i,j));
    end
    fprintf(f,'\r\n');
end

fclose(f);




































% e=210e9;
% g=84e9;
% a=2e-2;
% iy=1e-4;
% iz=2e-4;
% j=5e-5;
% eles=8;
% nodes=8;
% node_number=[1,2,3,4,5,6,7,8;
%     5,6,7,8,6,7,8,5];  %ÿ����Ԫ��Ӧ�ڵ�ı�š�
% node_coordinate=[0,0,0,0,5,0;
%                  0,0,4,0,5,4;
%                  4,0,4,4,5,4;
%                  4,0,0,4,5,0;
%                  0,5,0,0,5,4;
%                  0,5,4,4,5,4;
%                  4,5,4,4,5,0;
%                  4,5,0,0,5,0]';  %ÿ����Ԫ��Ӧ�ڵ�����ꡣ
% ele_info=zeros(10,eles);    %��Ԫ��Ϣ����ÿ��������e\g\a\iy\iz\j\l\cXx\cYx\cZx��
% ele_info(1,:)=e*ones(1,eles);
% ele_info(2,:)=a*ones(1,eles);
% ele_info(3,:)=iy*ones(1,eles);
% ele_info(4,:)=iz*ones(1,eles);
% ele_info(5,:)=j*ones(1,eles);
% ele_info(10,:)=g*ones(1,eles);
% ele_info(6:9,:)=frame3d_ele_info(node_coordinate);  %��ȫ��Ԫ��Ϣ�������Ϣ��
% r=zeros(12*eles,12);   %rΪ��Ԫת�ƾ���ļ��ϣ�ʮ��eles��12�С�
% for i=1:1:eles
%     r(12*i-11:12*i,:)=frame3d_r_mat(ele_info(7:9,i));
% end     %��ɵ�Ԫת�ƾ���ļ��ϡ�
% k_ele=zeros(12*eles,12);  %k_ele�ǵ�Ԫ�ھֲ�����ϵ�µĸնȾ���ʮ��eles��12�С�
% for i=1:1:eles
%     k_ele(12*i-11:12*i,:)=frame3d_k_ele_mat(ele_info(:,i));
% end     %��ɾֲ�����ϵ�µĵ�Ԫ�նȾ��󼯺ϡ�
% kk=zeros(12*eles,12); %kk�����������µ�Ԫ�նȾ���ļ��ϣ�12eles��12�С�
% for i=1:1:eles
%     kk(12*i-11:12*i,:)=frame3d_ele_stiff_mat(r(12*i-11:12*i,:),k_ele(12*i-11:12*i,:));
% end     %�����������ϵ�µ�Ԫ�նȾ��󼯺ϡ�
% k=zeros(nodes*6,nodes*6);   %k������նȾ���
% for i=1:1:eles
%     k=frame3d_assembly_mat(k,kk(12*i-11:12*i,:),node_number(:,i));
% end     %w�������նȾ���
%  m=k(25:48,25:48);
%  f=[zeros(12,1);-15e3;zeros(11,1)];
%  u=m\f;
%  U=[zeros(24,1);u]   %UΪ�ڵ���λ������
%  F=k*U   %FΪ�ڵ�����������
%  u=zeros(12,eles);   %uΪÿ��Ԫ���ڵ��λ�����С�
%  for i=1:1:eles
%     u(:,i)=[U(node_number(1,i)*6-5:node_number(1,i)*6);U(node_number(2,i)*6-5:node_number(2,i)*6)];
%  end
% f=zeros(12,eles);    %fΪÿ��Ԫ���ڵ���������У������¶�ÿ����Ԫ�����Ĳ���󣬻�Ӧ��ȥ��Ч�ڵ��غɡ�
% for i=1:1:eles
%     f(:,i)=k_ele(12*i-11:12*i,:)*r(12*i-11:12*i,:)*u(:,i);
% end
% f
