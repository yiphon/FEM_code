function y = tetra_solve_unknown( k,unknown_u_index,f,u,nodes )
%TETRA_SOLVE_UNKNOWN Summary of this function goes here
%   Detailed explanation goes here
%   利用分块矩阵计算未知的变量，力和位移。
answer={};
kk=k;
a=unknown_u_index;
n=size(a,1);
b=zeros(3*nodes-n,1);
j=1;
for i=1:1:3*nodes
    if ismember(i,a)==0
        b(j)=i;
        j=j+1;
    end
end
for i=1:1:n
    kk(i,:)=k(a(i),:);
end
for i=1:1:size(b,1)
    kk(n+i,:)=k(b(i),:);
end
kkk=kk;
for i=1:1:n
    kk(:,i)=kkk(:,a(i));
end
for i=1:1:size(b,1)
    kk(:,n+i)=kkk(:,b(i));
end     %需要先后进行矩阵的行变换和列变换。
k1=kk(1:n,1:n);
k2=kk(1:n,n+1:3*nodes);
k3=kk(n+1:3*nodes,1:n);
k4=kk(n+1:3*nodes,n+1:3*nodes);
u_solved=k1\(f-k2*u);
f_solved=k3*u_solved+k4*u;
answer{1,1}=u_solved;
answer{1,2}=f_solved;
F=zeros(3*nodes,1);
U=zeros(3*nodes,1);
index_f_unknown=1;
index_u_known=1;
index_f_known=1;
index_u_unknown=1;
for i=1:1:3*nodes
    if ismember(i,a)==0     %位移为已知量，力为未知量的编号。
        F(i)=f_solved(index_f_unknown);
        index_f_unknown=index_f_unknown+1;
        U(i)=u(index_u_known);
        index_u_known=index_u_known+1;
    else    %位移为未知量，力为已知量的编号。
        F(i)=f(index_f_known);
        index_f_known=index_f_known+1;
        U(i)=u_solved(index_u_unknown);
        index_u_unknown=index_u_unknown+1;
    end
end
answer{1,3}=U;
answer{1,4}=F;
answer{1,5}=k1;
answer{1,6}=k2;
answer{1,7}=k3;
answer{1,8}=k4;
y=answer;


end

