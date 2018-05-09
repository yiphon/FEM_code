%立体刚架元算例。
e=207e9;
g=84e9;
a=18e-4;
iy=0.06*0.03^3/12;
iz=0.03*0.06^3/12;
j=5e-5;
eles=1;
nodes=2;
node_number=[1;2];  %每个单元对应节点的编号。
node_coordinate=[0,0,0,1.5,0,0]';  %每个单元对应节点的坐标。
ele_info=zeros(10,eles);    %单元信息矩阵，每列依次是e\a\iy\iz\j\l\cXx\cYx\cZx\g。
ele_info(1,:)=e*ones(1,eles);
ele_info(2,:)=a*ones(1,eles);
ele_info(3,:)=iy*ones(1,eles);
ele_info(4,:)=iz*ones(1,eles);
ele_info(5,:)=j*ones(1,eles);
ele_info(10,:)=g*ones(1,eles);
ele_info(6:9,:)=frame3d_ele_info(node_coordinate);  %补全单元信息矩阵的信息。
r=zeros(12*eles,12);   %r为单元转移矩阵的集合，十二eles行12列。
for i=1:1:eles
    r(12*i-11:12*i,:)=frame3d_r_mat(ele_info(7:9,i));
end     %完成单元转移矩阵的集合。
k_ele=zeros(12*eles,12);  %k_ele是单元在局部坐标系下的刚度矩阵，十二eles行12列。
for i=1:1:eles
    k_ele(12*i-11:12*i,:)=frame3d_k_ele_mat(ele_info(:,i));
end     %完成局部坐标系下的单元刚度矩阵集合。
kk=zeros(12*eles,12); %kk是整体坐标下单元刚度矩阵的集合，12eles行12列。
for i=1:1:eles
    kk(12*i-11:12*i,:)=frame3d_ele_stiff_mat(r(12*i-11:12*i,:),k_ele(12*i-11:12*i,:));
end     %完成整体坐标系下单元刚度矩阵集合。
k=zeros(nodes*6,nodes*6);   %k是总体刚度矩阵。
for i=1:1:eles
    k=frame3d_assembly_mat(k,kk(12*i-11:12*i,:),node_number(:,i));
end     %w完成总体刚度矩阵。
 m=k(7:12,7:12);
 f=[0;500;300;0;0;0];
 u=m\f;
 U=[zeros(6,1);u]   %U为节点总位移列阵。
 F=k*U   %F为节点总受力列阵。
 u=zeros(12,eles);   %u为每单元两节点的位移阵列。
 for i=1:1:eles
    u(:,i)=[U(node_number(1,i)*6-5:node_number(1,i)*6);U(node_number(2,i)*6-5:node_number(2,i)*6)];
 end
f=zeros(12,eles);    %f为每单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
for i=1:1:eles
    f(:,i)=k_ele(12*i-11:12*i,:)*r(12*i-11:12*i,:)*u(:,i);
end
f
max_bending_sigma=sym(zeros(2,eles));   %用以求各单元的最大弯曲应力值。
for i=1:1:eles
    max_bending_sigma(:,i)=frame3d_max_bending_sigma(ele_info(1,i),ele_info(6,i),u(:,i));
end
max_bending_sigma=vpa(max_bending_sigma,5)
