%立体刚架元算例。

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
node_number=s(8:7+eles,1:2)';   %每个单元对应节点的编号。
node_coordinate=s(8+eles:7+eles+nodes,1:3)';
node_coordinate=truss3d_ele_coordinate(node_coordinate,node_number);   %每个单元对应节点的坐标。
unknown_u_index=s(8+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';
f_known=s(9+eles+nodes,1:size(unknown_u_index,1))';
u_known=s(10+eles+nodes,1:6*nodes-size(unknown_u_index,1))';

ele_info=zeros(10,eles);    %单元信息矩阵，每列依次是e\a\iy\iz\j\l\cXx\cYx\cZx\g。
ele_info(1,:)=e;
ele_info(2,:)=a;
ele_info(3,:)=iy;
ele_info(4,:)=iz;
ele_info(5,:)=j;
ele_info(10,:)=g;
ele_info(6:9,:)=frame3d_ele_info(node_coordinate);  %补全单元信息矩阵的信息。

r=zeros(12,12,eles);   %r为单元转移矩阵的集合，每单元十二行12列。
for i=1:1:eles
    r(:,:,i)=frame3d_r_mat(ele_info(7:9,i));
end     %完成单元转移矩阵的集合。

k_ele=zeros(12,12,eles);  %k_ele是单元在局部坐标系下的刚度矩阵，每单元十二行12列。
for i=1:1:eles
    k_ele(:,:,i)=frame3d_k_ele_mat(ele_info(:,i));
end     %完成局部坐标系下的单元刚度矩阵集合。

kk=zeros(12,12,eles); %kk是整体坐标下单元刚度矩阵的集合，每单元12行12列。
for i=1:1:eles
    kk(:,:,i)=frame3d_ele_stiff_mat(r(:,:,i),k_ele(:,:,i));
end     %完成整体坐标系下单元刚度矩阵集合。


k=zeros(nodes*6,nodes*6);   %k是总体刚度矩阵。
for i=1:1:eles
    k=frame3d_assembly_mat(k,kk(:,:,i),node_number(:,i));
end     %w完成总体刚度矩阵。

unknown_var=frame3d_solve_unknown(k,unknown_u_index,f_known,u_known,nodes);     %利用分块矩阵计算未知的变量，力和位移。
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

u=zeros(12,eles);   %u为每单元两节点的位移阵列。
for i=1:1:eles
   u(:,i)=[U(node_number(1,i)*6-5:node_number(1,i)*6);U(node_number(2,i)*6-5:node_number(2,i)*6)];
end

f_nodal=zeros(12,eles);    %f为每单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
for i=1:1:eles
    f_nodal(:,i)=k_ele(:,:,i)*r(:,:,i)*u(:,i);
end


%以下写文件。
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
%     5,6,7,8,6,7,8,5];  %每个单元对应节点的编号。
% node_coordinate=[0,0,0,0,5,0;
%                  0,0,4,0,5,4;
%                  4,0,4,4,5,4;
%                  4,0,0,4,5,0;
%                  0,5,0,0,5,4;
%                  0,5,4,4,5,4;
%                  4,5,4,4,5,0;
%                  4,5,0,0,5,0]';  %每个单元对应节点的坐标。
% ele_info=zeros(10,eles);    %单元信息矩阵，每列依次是e\g\a\iy\iz\j\l\cXx\cYx\cZx。
% ele_info(1,:)=e*ones(1,eles);
% ele_info(2,:)=a*ones(1,eles);
% ele_info(3,:)=iy*ones(1,eles);
% ele_info(4,:)=iz*ones(1,eles);
% ele_info(5,:)=j*ones(1,eles);
% ele_info(10,:)=g*ones(1,eles);
% ele_info(6:9,:)=frame3d_ele_info(node_coordinate);  %补全单元信息矩阵的信息。
% r=zeros(12*eles,12);   %r为单元转移矩阵的集合，十二eles行12列。
% for i=1:1:eles
%     r(12*i-11:12*i,:)=frame3d_r_mat(ele_info(7:9,i));
% end     %完成单元转移矩阵的集合。
% k_ele=zeros(12*eles,12);  %k_ele是单元在局部坐标系下的刚度矩阵，十二eles行12列。
% for i=1:1:eles
%     k_ele(12*i-11:12*i,:)=frame3d_k_ele_mat(ele_info(:,i));
% end     %完成局部坐标系下的单元刚度矩阵集合。
% kk=zeros(12*eles,12); %kk是整体坐标下单元刚度矩阵的集合，12eles行12列。
% for i=1:1:eles
%     kk(12*i-11:12*i,:)=frame3d_ele_stiff_mat(r(12*i-11:12*i,:),k_ele(12*i-11:12*i,:));
% end     %完成整体坐标系下单元刚度矩阵集合。
% k=zeros(nodes*6,nodes*6);   %k是总体刚度矩阵。
% for i=1:1:eles
%     k=frame3d_assembly_mat(k,kk(12*i-11:12*i,:),node_number(:,i));
% end     %w完成总体刚度矩阵。
%  m=k(25:48,25:48);
%  f=[zeros(12,1);-15e3;zeros(11,1)];
%  u=m\f;
%  U=[zeros(24,1);u]   %U为节点总位移列阵。
%  F=k*U   %F为节点总受力列阵。
%  u=zeros(12,eles);   %u为每单元两节点的位移阵列。
%  for i=1:1:eles
%     u(:,i)=[U(node_number(1,i)*6-5:node_number(1,i)*6);U(node_number(2,i)*6-5:node_number(2,i)*6)];
%  end
% f=zeros(12,eles);    %f为每单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
% for i=1:1:eles
%     f(:,i)=k_ele(12*i-11:12*i,:)*r(12*i-11:12*i,:)*u(:,i);
% end
% f
