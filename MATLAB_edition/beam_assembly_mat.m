function y = beam_assembly_mat( K,kk,node )
%BEAM_ASSEMBLY_MAT Summary of this function goes here
%   Detailed explanation goes here
%   �����ܸնȾ���kkΪ��λ�նȾ���node_numberΪ[i;j]��ʽ�Ľڵ��š�
i=node(1);
j=node(2);
a=[2*i-1;2*i;2*j-1;2*j];
for m=1:1:4
    for n=1:1:4
        K(a(m),a(n))=K(a(m),a(n))+kk(m,n);
    end
end
y=K;
end

