%平面刚架元算例。

s=xlsread('frame2d_eg_input.xlsx');
s(all(isnan(s),2),:)=[];
eles=s(1,1);
nodes=s(1,2);
e=s(2,1:eles);
a=s(3,1:eles);
iz=s(4,1:eles);
node_number=s(5:4+eles,1:2)';   %每个单元对应节点的编号。
node_coordinate=s(5+eles:4+eles+nodes,1:2)';
node_coordinate=truss2d_ele_coordinate(node_coordinate,node_number);   %每个单元对应节点的坐标。
unknown_u_index=s(5+eles+nodes,:);
unknown_u_index(:,all(isnan(unknown_u_index),1))=[];
unknown_u_index=unknown_u_index';
f_nodal=s(6+eles+nodes,1:size(unknown_u_index,1))';
u_known=s(7+eles+nodes,1:3*nodes-size(unknown_u_index,1))';
dis_f_val=s(8+eles+nodes,1);
dis_f_ele=s(8+eles+nodes,2);
dis_f_range=s(8+eles+nodes,3);

ele_info=zeros(6,eles);    %单元信息矩阵，每列依次是e\a\iz\l\cos\sin。
ele_info(1,:)=e;
ele_info(2,:)=a;
ele_info(3,:)=iz;
ele_info(4:6,:)=frame2d_ele_info(node_coordinate);  %补全单元信息矩阵的信息。

r=zeros(6,6,eles);   %r为单元转移矩阵的集合，每单元六行六列。
for i=1:1:eles
    r(:,:,i)=frame2d_r_mat(ele_info(5:6,i));
end     %完成单元转移矩阵的集合。


k_ele=zeros(6,6,eles);  %k_ele单元在局部坐标系下的刚度矩阵，每单元六行六列。
for i=1:1:eles
    k_ele(:,:,i)=frame2d_k_ele_mat(ele_info(1:4,i));
end     %完成局部坐标系下的单元刚度矩阵集合。
k_ele
kk=zeros(6,6,eles); %kk是整体坐标下单元刚度矩阵的集合，每单元六行六列。
for i=1:1:eles
    kk(:,:,i)=frame2d_ele_stiff_mat(r(:,:,i),k_ele(:,:,i));
end     %完成整体坐标系下单元刚度矩阵集合。


k=zeros(nodes*3,nodes*3);   %k是总体刚度矩阵。
for i=1:1:eles
    k=frame2d_assembly_mat(k,kk(:,:,i),node_number(:,i));
end     %w完成总体刚度矩阵。


equi=beam_work_equi_distributed(dis_f_range,dis_f_val);   %equi为等效节点载荷列阵。
f_known=f_nodal;
f_known(find(unknown_u_index==3*node_number(1,dis_f_ele)-1))= ...
f_known(find(unknown_u_index==3*node_number(1,dis_f_ele)-1))+equi(1);
f_known(find(unknown_u_index==3*node_number(1,dis_f_ele)))= ...
f_known(find(unknown_u_index==3*node_number(1,dis_f_ele)))+equi(2);
f_known(find(unknown_u_index==3*node_number(2,dis_f_ele)-1))= ...
f_known(find(unknown_u_index==3*node_number(2,dis_f_ele)-1))+equi(3);
f_known(find(unknown_u_index==3*node_number(2,dis_f_ele)))= ...
f_known(find(unknown_u_index==3*node_number(2,dis_f_ele)))+equi(4);

unknown_var=tetra_solve_unknown(k,unknown_u_index,f_known,u_known,nodes);     %利用分块矩阵计算未知的变量，力和位移。
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

u=zeros(6,eles);   %u为每单元两节点的位移阵列。
for i=1:1:eles
    u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
end
f_nodal=zeros(6,eles);    %f为每单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
for i=1:1:eles
    f_nodal(:,i)=k_ele(:,:,i)*r(:,:,i)*u(:,i);
    if i==dis_f_ele
        equi=[0;equi(1);equi(2);0;equi(3);equi(4)];
        f_nodal(:,i)=f_nodal(:,i)-equi;
    end
end

for i=1:1:eles
    beam_deflection_plot(ele_info(4,i),u([2,3,5,6],i),eles,i);
end     %绘出各单元挠度转角图像。


%以下写文件。
f=fopen('frame2d_eg_output.xls','w');
   
fprintf(f, 'element stiffness matrices:\r\n\r\n');   
for i=1:eles
    for j=1:6
       for m=1:6
           fprintf(f,'%8.4f\t',kk(j,m,i));
       end
       fprintf(f, '\r\n');    
    end
    fprintf(f,'\r\n')
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

fprintf(f,'\r\n\r\nf_nodal= \r\n\r\n');
for i=1:size(f_nodal,1)
    for j=1:size(f_nodal,2)
        fprintf(f,'%8.4f\t\t',f_nodal(i,j));
    end
    fprintf(f,'\r\n');
end

fclose(f);





%  e=210e9;
%  a=1e-2;
%  iz=9e-5;
%  eles=3;
%  nodes=4;
%  node_number=[1,2,3;2,3,4];  %每个单元对应节点的编号。
%  node_coordinate=[0,2,7;
%      0,3,3;
%      2,7,9;
%      3,3,0];  %每个单元对应节点的坐标。
%  ele_info=zeros(6,eles);    %单元信息矩阵，每列依次是e\a\iz\l\cos\sin。
%  ele_info(1,:)=e*ones(1,eles);
%  ele_info(2,:)=a*ones(1,eles);
%  ele_info(3,:)=iz*ones(1,eles);
%  ele_info(4:6,:)=frame2d_ele_info(node_coordinate);  %补全单元信息矩阵的信息。
% r=zeros(6*eles,6);   %r为单元转移矩阵的集合，六eles行六列。
%  for i=1:1:eles
%      r(6*i-5:6*i,:)=frame2d_r_mat(ele_info(5:6,i));
%  end     %完成单元转移矩阵的集合。
%  k_ele=zeros(6*eles,6);  %k_ele是单元在局部坐标系下的刚度矩阵，六eles行六列。
%  for i=1:1:eles
%      k_ele(6*i-5:6*i,:)=frame2d_k_ele_mat(ele_info(1:4,i));
%  end     %完成局部坐标系下的单元刚度矩阵集合。
%  kk=zeros(6*eles,6); %kk是整体坐标下单元刚度矩阵的集合，六eles行六列。
%  for i=1:1:eles
%      kk(6*i-5:6*i,:)=frame2d_ele_stiff_mat(r(6*i-5:6*i,:),k_ele(6*i-5:6*i,:));
%  end     %完成整体坐标系下单元刚度矩阵集合。
%  k=zeros(nodes*3,nodes*3);   %k是总体刚度矩阵。
%  for i=1:1:eles
%      k=frame2d_assembly_mat(k,kk(6*i-5:6*i,:),node_number(:,i));
%  end     %w完成总体刚度矩阵。
%  k
%  
%   m=k(4:9,4:9);
%   equi=beam_work_equi_distributed(5,-5000);   %equi为等效节点载荷列阵。
%   f=[20000;equi(1);equi(2);0;equi(3);equi(4)];
%   u=m\f;
%   U=[0;0;0;u;0;0;0]   %U为节点总位移列阵。
%   F=k*U   %F为节点总受力列阵。
%   
%   u=zeros(6,eles);   %u为每单元两节点的位移阵列。
%   for i=1:1:eles
%      u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
%   end
%  f=zeros(6,eles);    %f为每单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
%  for i=1:1:eles
%      f(:,i)=k_ele(6*i-5:6*i,:)*r(6*i-5:6*i,:)*u(:,i);
%      if i==2
%          equi=[0;equi(1);equi(2);0;equi(3);equi(4)];
%          f(:,i)=f(:,i)-equi;
%      end
%  end
%  f
%  
%  for i=1:1:eles
%      beam_deflection_plot(ele_info(4,i),u([2,3,5,6],i),eles,i);
%  end     %绘出各单元挠度转角图像。

 
 



