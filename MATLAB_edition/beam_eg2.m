%  梁元算法实例。
e=70e9;
I=40e-6;
ele_info=[e,e,;I,I;3,3];   %生成单元信息矩阵。
eles_beam=2;
nodes_beam=3;
node_number=[1,2;2,3];
kk=zeros(4*eles_beam,4); %kk是各单元刚度矩阵的集合。
for i=1:1:eles_beam
    kk(4*i-3:4*i,:)=beam_ele_stiff_mat(ele_info(:,i));
end     %生成单元刚度矩阵。
kk_bar=5e6*[1 -1;-1 1]; %kk_bar是弹簧的刚度矩阵。
k=zeros(nodes_beam*2+1,nodes_beam*2+1);   %k为整体刚度矩阵。
for i=1:1:eles_beam
    k=beam_assembly_mat(k,kk(4*i-3:4*i,:),node_number(:,i));
end
k(3,3)=k(3,3)+kk_bar(1,1);
k(3,7)=k(3,7)+kk_bar(1,2);
k(7,3)=k(7,3)+kk_bar(2,1);
k(7,7)=k(7,7)+kk_bar(2,2);  %生成整体刚度矩阵。
% equi=zeros(4,2);    %equi为等效节点载荷矩阵。
% equi(:,1)=beam_work_equi_distributed(3);
% equi(:,2)=beam_work_equi_central(2,4,-30e3);
% f=[equi(4,1);equi(2,2);equi(4,2)];
m=k([3,4,6],[3,4,6]);
f=[-1e4;0;0];
u=m\f;  %求出位置的三个位移。
U=[0;0;u(1);u(2);0;u(3);0]   %U为总位移矩阵。
F=k*U   %F为总节点力矩阵。
u=zeros(4,eles_beam);    %u为各单元两节点的位移矩阵。
for i=1:1:eles_beam
    u(:,i)=[U(node_number(1,i)*2-1:node_number(1,i)*2);U(node_number(2,i)*2-1:node_number(2,i)*2)];
end
for i=1:1:eles_beam
   beam_deflection_plot(ele_info(3,i),u(:,i),3,i);
end     %绘出各单元挠度与转角图像。
f_beam=zeros(4,eles_beam);    %f为各单元两节点的力矩阵。
for i=1:1:eles_beam
    f_beam(:,i)=beam_node_force(kk(4*i-3:4*i,:),u(:,i));
end
f_bar=kk_bar*[u(3);u(7)];
f_beam
f_bar



    