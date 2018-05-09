%平面刚架元算例。
e=70e9;
a=1e-2;
iz=1e-5;
eles_beam=1;
nodes_beam=2;
node_number=[1;2];  %每个单元对应节点的编号。
node_coordinate=[0;0;4;0];  %每个单元对应节点的坐标。
ele_info=zeros(6,eles_beam);    %单元信息矩阵，每列依次是e\a\iz\l\cos\sin。
ele_info(1,:)=e*ones(1,eles_beam);
ele_info(2,:)=a*ones(1,eles_beam);
ele_info(3,:)=iz*ones(1,eles_beam);
ele_info(4:6,:)=frame2d_ele_info(node_coordinate);  %补全单元信息矩阵的信息。
r=zeros(6*eles_beam,6);   %r为单元转移矩阵的集合，六eles行六列。
for i=1:1:eles_beam
    r(6*i-5:6*i,:)=frame2d_r_mat(ele_info(5:6,i));
end     %完成单元转移矩阵的集合。
k_ele=zeros(6*eles_beam,6);  %k_ele是单元在局部坐标系下的刚度矩阵，六eles行六列。
for i=1:1:eles_beam
    k_ele(6*i-5:6*i,:)=frame2d_k_ele_mat(ele_info(1:4,i));
end     %完成局部坐标系下的单元刚度矩阵集合。
kk=zeros(6*eles_beam,6); %kk是整体坐标下单元刚度矩阵的集合，六eles行六列。
for i=1:1:eles_beam
    kk(6*i-5:6*i,:)=frame2d_ele_stiff_mat(r(6*i-5:6*i,:),k_ele(6*i-5:6*i,:));
end     %完成整体坐标系下单元刚度矩阵集合。
k=zeros(nodes_beam*3+2,nodes_beam*3+2);   %k是总体刚度矩阵。
for i=1:1:eles_beam
    k=frame2d_assembly_mat(k,kk(6*i-5:6*i,:),node_number(:,i));
end     %w完成总体刚度矩阵刚架部分。
k_truss=ele_stiffness_mat_2dtruss(1,1,0,0,4,3);
k_truss=5*5e6*k_truss;
k([1,2,7,8],[1,2,7,8])=k([1,2,7,8],[1,2,7,8])+k_truss;  %完成总体刚度矩阵桁架部分。

 m=k(1:3,1:3);
%  equi=beam_work_equi_distributed(5,-5000);   %equi为等效节点载荷列阵。
 %f=[20000;equi(1);equi(2);0;equi(3);equi(4)];
 f=[0;-1e4;0];
 u=m\f;
 U=[u;0;0;0;0;0]   %U为节点总位移列阵。
 F=k*U   %F为节点总受力列阵。
 
 u_frame=zeros(6,eles_beam);   %u为每刚架单元两节点的位移阵列。
 for i=1:1:eles_beam
    u_frame(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
 end
f_frame=zeros(6,eles_beam);    %f为每刚架单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
for i=1:1:eles_beam
    f_frame(:,i)=k_ele(6*i-5:6*i,:)*r(6*i-5:6*i,:)*u_frame(:,i);    
end
f_frame

u_truss=U([1,2,7,8]);
f_truss=ele_inner_force_2dtruss(1,1,0,0,4,3,u_truss);
f_truss=5*5e6*f_truss;
f_truss

for i=1:1:eles_beam
    beam_deflection_plot(ele_info(4,i),u_frame([2,3,5,6],i),eles_beam,i);
end     %绘出各单元挠度转角图像。
