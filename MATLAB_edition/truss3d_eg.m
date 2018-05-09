%立体桁架结构有限元算法实例。
eles_str=['Input the number of the truss elements:' 10];
eles=input(eles_str)     %eles表示结构的单元个数。
nodes_str=['Input the number of nodes in the structure:' 10];
nodes=input(nodes_str)      %nodes表示结构中的节点个数

e_str=['Input the elastic modulus for each elements, split with comma:' 10];
E=input(e_str,'s');
E=str2num(E)        %E为各单元弹性模量的集合。
a_str=['Input the area of transverse cross section for each elements, split with comma:' 10];
A=input(a_str,'s');
A=str2num(A)        %A为各单元横截面积的集合。

n_n_str=['Input the node sequences for each elements, split with comma, ' 10 'and ' ...
            'separate the elements with semicolon. ' 10 'eg:  2,3;4,5 then press Enter.' 10];
node_number=input(n_n_str,'s');
node_number=str2num(node_number)'   %node_number为各单元对应节点的集合，每列对应一个单元。

n_c_str=['Input the coordinate of each nodes in sequence, split with comma, ' 10 'and ' ...
            'separate the nodes with semicolon. ' 10 'eg:  0,0,0;0.5,2.5,3.2 then press Enter.' 10];
node_coordinate=input(n_c_str,'s');
node_coordinate=str2num(node_coordinate)'       %node_coordinate为各节点的坐标集合，每列对应一个一个节点。
node_coordinate=truss3d_ele_coordinate(node_coordinate,node_number);  %node_coordinate生成每个单元对应节点的坐标列阵集合。

ele_info=zeros(4,eles);
 for i=1:1:eles
     ele_info(:,i)=truss3d_ele_info(node_coordinate(:,i));
 end     %完成单元基本信息矩阵，包含方向余弦与长度。

 kk=zeros(6,6,eles);   %单元刚度矩阵的集合。
 for i=1:1:eles
     kk(:,:,i)=truss3d_ele_stiffness_mat(E(i),A(i),ele_info(:,i));
 end     %完成单元刚度矩阵。
 kk
 
 k=zeros(3*nodes,3*nodes);
 for i=1:1:eles
     k=truss3d_ele_mat_assemble(k,kk(:,:,i),node_number(:,i));
 end     %装配整体刚度矩阵。
 k
 
u_u_i_str=['Input the indexs of the unknown nodal displacement in the U vector ' 10 'of the ' ...
    'global stiffness equations K*U=F, split with semicolon.' 10]; 
unknown_u_index=input(u_u_i_str,'s');
unknown_u_index=str2num(unknown_u_index)       %unknown_u_index是总体刚度矩阵中未知节点位移序号的索引。

k_u_str=['Input the known nodal displacement sequence in the U vector ' 10 'of the ' ...
    'global stiffness equations K*U=F, split with semicolon.' 10]; 
known_u=input(k_u_str,'s');
known_u=str2num(known_u)        %known_u是总体刚度方程中已知节点位移的集合。


k_f_str=['Input the known nodal force sequence in the F vector ' 10 'of the ' ...
    'global stiffness equations K*U=F, split with semicolon.' 10]; 
known_f=input(k_f_str,'s');
known_f=str2num(known_f)       %known_f是总体刚度方程中已知节点力的集合。

unknown_var=tetra_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
unknown_u=unknown_var{1,1}
unknown_f=unknown_var{1,2}
U=unknown_var{1,3}  %U为总位移列阵。
F=unknown_var{1,4}   %F为总受力列阵。

u=zeros(6,eles);    %u为各单元对应节点的位移列阵集合。
 for i=1:1:eles
       u(:,i)=[U(node_number(1,i)*3-2:node_number(1,i)*3);U(node_number(2,i)*3-2:node_number(2,i)*3)];
 end

 force=zeros(eles,1);
 stress=zeros(eles,1);
 for i=1:1:eles
     force(i)=truss3d_ele_inner_force(E(i),A(i),ele_info(:,i),u(:,i));
     stress(i)=force(i)/A(i);
 end
 force
 stress       %求各单元的内力与应力。


 

%   E= 2e11;
%   A=[0.003;0.003;0.003;0.003];
%   eles=4;
%   nodes=5
%   node_number=[1,2,3,4;5,5,5,5];  %每个单元对应起点与终点的节点编号列阵。
%   node_coordinate=[0,0,-3,0,5,0;
%                    -3,0,0,0,5,0;
%                   0,0,3,0,5,0;
%                    4,0,0,0,5,0]';       %每行表示每个单元起点与终点的节点坐标序列。
%   ele_info=zeros(4,eles);
%   for i=1:1:eles
%      ele_info(:,i)=truss3d_ele_info(node_coordinate(:,i));
%   end     %完成单元基本信息矩阵，包含方向余弦与长度。
%  kk=zeros(6*eles,6);   %单元刚度矩阵的集合。
%   for i=1:1:eles
%       kk(6*i-5:6*i,:)=truss3d_ele_stiffness_mat(E,A(i),ele_info(:,i));
%   end     %完成单元刚度矩阵。
%   k=zeros(3*nodes,3*nodes);
%   for i=1:1:eles
%       k=truss3d_ele_mat_assemble(k,kk(6*i-5:6*i,:),node_number(:,i));
%   end     %装配整体刚度矩阵。
%   m=k(13:15,13:15);
%   f=[15e3;0;-2e4];
%   u=m\f;          %根据已知条件删减整体刚度矩阵，以便求得最终所有的未知量。
%   U=[0;0;0;0;0;0;0;0;0;0;0;0;u]        %U为最终位移矩阵。
%   F=k*U      %F为最终力矩阵。
%   u=zeros(6,eles);
%   for i=1:1:eles
%      u(:,i)=[U((3*node_number(1,i)-2):(3*node_number(1,i)));U((3*node_number(2,i)-2):(3*node_number(2,i)))];
%   end         %每个单元对应的两个节点的位移分量。
%   force=zeros(eles,1);
%   stress=zeros(eles,1);
%   for i=1:1:eles
%       force(i)=truss3d_ele_inner_force(E,A(i),ele_info(:,i),u(:,i));
%       stress(i)=force(i)/A(i);
%   end
%   force
%   stress       %求各单元的内力与应力。
% 
% 
% 
% 
% 
% 
% 
