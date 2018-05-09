%平面桁架结构有限元算法实例。
eles_str=['Input the number of the plane truss elements:' 10];
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
            'separate the nodes with semicolon. ' 10 'eg:  0,0;0.5,2.5 then press Enter.' 10];
node_coordinate=input(n_c_str,'s');
node_coordinate=str2num(node_coordinate)'       %node_coordinate为各节点的坐标集合，每列对应一个一个节点。
node_coordinate=truss2d_ele_coordinate(node_coordinate,node_number);  %node_coordinate生成每个单元对应节点的坐标列阵集合。

kk=zeros(4,4,eles);     %kk为各单元刚度矩阵的集合。
for i=1:1:eles
    kk(:,:,i)=truss2d_ele_stiff_mat(E(i),A(i),node_coordinate(:,i));
end
kk

k=zeros(2*nodes,2*nodes);   %k为总体的刚度矩阵。
for i=1:1:eles
    k=truss2d_ele_mat_assembly(k,kk(:,:,i),node_number(:,i));
end
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

unknown_var=tri_solve_unknown(k,unknown_u_index,known_f,known_u,nodes);     %利用分块矩阵计算未知的变量，力和位移。
unknown_u=unknown_var{1,1}
unknown_f=unknown_var{1,2}
U=unknown_var{1,3}  %U为总位移列阵。
F=unknown_var{1,4}   %F为总受力列阵。

u=zeros(4,eles);    %u为各单元对应节点的位移列阵集合。
 for i=1:1:eles
       u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2)];
 end

 force=zeros(eles,1);   %force为各单元所受的内力,stress对应应力。
 stress=zeros(eles,1);
 for i=1:1:eles
     force(i)=truss2d_ele_inner_force(E(i),A(i),u(:,i),node_coordinate(:,i));
     stress(i)=force(i)/A(i);
 end
 force
 stress













% E= 2.95e11;
% A=0.0001;
% a=ele_stiffness_mat_2dtruss(E,A,0,0,0.4,0);
% b=ele_stiffness_mat_2dtruss(E,A,0.4,0.3,0.4,0);
% c=ele_stiffness_mat_2dtruss(E,A,0,0,0.4,0.3);
% d=ele_stiffness_mat_2dtruss(E,A,0,0.3,0.4,0.3);   %产生单元刚度矩阵。
% k=zeros(8);
% k=ele_mat_assemble_2dtruss(k,a,1,2);
% k=ele_mat_assemble_2dtruss(k,b,3,2);
% k=ele_mat_assemble_2dtruss(k,c,1,3);
% k=ele_mat_assemble_2dtruss(k,d,4,3);              %装配整体刚度矩阵。
% m=[k(3,3),k(3,5:6);k(5:6,3),k(5:6,5:6)];          %根据已知条件删减整体刚度矩阵，以便求得最终所有的未知量。
% f=[20e3;2.8e3;-25e3];
% u=m\f;
% U=[0;0;u(1);0;u(2);u(3);0;0.0001]  %U为最终位移矩阵。
% F=k*U                         %F为最终力矩阵。
% u1=U(1:4);
% force1=ele_inner_force_2dtruss(E,A,0,0,0.4,0,u1)    %求各单元的内力与应力。
% stress1=ele_inner_stress_2dtruss(E,0,0,0.4,0,u1)
% u2=[U(5);U(6);U(3);U(4)];
% force2=ele_inner_force_2dtruss(E,A,0.4,0.3,0.4,0,u2)
% stress2=ele_inner_stress_2dtruss(E,0.4,0.3,0.4,0,u2)
% u3=[U(1);U(2);U(5);U(6)];
% force3=ele_inner_force_2dtruss(E,A,0,0,0.4,0.3,u3)
% u4=[U(7);U(8);U(5);U(6)];
% force4=ele_inner_force_2dtruss(E,A,0,0.3,0.4,0.3,u4)