%平面刚架元算例。
e=1e7;
a=8;
iz=2*4^3/12;
eles=2;
nodes=3;
node_number=[3,2;1,2]';  %每个单元对应节点的编号。
node_coordinate=[0,30,30,30;0,0,30,30]';  %每个单元对应节点的坐标。
ele_info=zeros(6,eles);    %单元信息矩阵，每列依次是e\a\iz\l\cos\sin。
ele_info(1,:)=e*ones(1,eles);
ele_info(2,:)=a*ones(1,eles);
ele_info(3,:)=iz*ones(1,eles);
ele_info(4:6,:)=frame2d_ele_info(node_coordinate);  %补全单元信息矩阵的信息。
r=zeros(6*eles,6);   %r为单元转移矩阵的集合，六eles行六列。
for i=1:1:eles
    r(6*i-5:6*i,:)=frame2d_r_mat(ele_info(5:6,i));
end     %完成单元转移矩阵的集合。
k_ele=zeros(6*eles,6);  %k_ele是单元在局部坐标系下的刚度矩阵，六eles行六列。
for i=1:1:eles
    k_ele(6*i-5:6*i,:)=frame2d_k_ele_mat(ele_info(1:4,i));
end     %完成局部坐标系下的单元刚度矩阵集合。
kk=zeros(6*eles,6); %kk是整体坐标下单元刚度矩阵的集合，六eles行六列。
for i=1:1:eles
    kk(6*i-5:6*i,:)=frame2d_ele_stiff_mat(r(6*i-5:6*i,:),k_ele(6*i-5:6*i,:));
end     %完成整体坐标系下单元刚度矩阵集合。
kk

k=zeros(nodes*3,nodes*3);   %k是总体刚度矩阵。
for i=eles:-1:1
    k=frame2d_assembly_mat(k,kk(6*i-5:6*i,:),node_number(:,i));
end     %w完成总体刚度矩阵。

k
equi=beam_work_equi_distributed(30,-20); 
equi
unknown_u_index=[4,5,6]';   %unknown_u_index是未知位移序号的索引。
 known_f=[0;-1200+equi(3);-1500+equi(4)];
 known_u=zeros(6,1);
 unknown_var=frame2d_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
 unknown_u=unknown_var{1,1}
 unknown_f=unknown_var{1,2}
U=unknown_var{1,3}  %U为总位移列阵。
F=unknown_var{1,4}   %F为总受力列阵。
F([8,9,5,6])=F([8,9,5,6])-equi;
f_node_reaction=F([5,6,8,9])


 u=zeros(6,eles);   %u为每单元两节点的位移阵列。
 for i=1:1:eles
    u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
 end
f=zeros(6,eles);    %f为每单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
for i=1:1:eles
    f(:,i)=k_ele(6*i-5:6*i,:)*r(6*i-5:6*i,:)*u(:,i);
    if i==1
        equi=[0;equi(1);equi(2);0;equi(3);equi(4)];
        f(:,i)=f(:,i)-equi;
    end
end
f


max_sigma=sym(zeros(2,eles));   %用以求各单元的最大应力值。
for i=1:1:eles
    max_sigma(:,i)=frame2d_max_sigma(ele_info(1,i),ele_info(4,i),u(:,i),f(:,i),ele_info(2,i));
end
max_sigma=vpa(max_sigma,5)


 
 



