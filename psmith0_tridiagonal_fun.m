function [x,xs,g]=psmith0_tridiagonal_fun(A,B)
% AOE4404HW2P1 by Patrick Smith (psmith0), 2/22/13
%
% This function solves for a set of linear equations
% when the A matrix is aTridiagonal using Thomas algorithm
%
%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%
% A: Tridiagonal A matrix
% B: B vector in Ax=B
%
%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%
% x: Solution vector x in Ax=B
% xs: Solution vector in L*xstar=B
% g: Elements in U matrix above diagonal

%% preallocate
otd=0;
[m,n]=size(A);
g=zeros(m,1);xs=g;x=g;

%% check for things that will break code:

if m~=n
    display('A matrix is not Square');
    x=[];    xs=[];    g=[];
    return;
end

if min(size(B))>1
    display('B matrix not singleton')
    x=[];    xs=[];    g=[];
    return;
end

if max(size(B))~=size(A,1)
    display('dimensions of A and B do not match')
    x=[];    xs=[];    g=[];
    return;
end

%% Check that the matrix is tridiagonal, solve for gamma, xstar

for i=1:size(A,1)
    if i<3;
       otd=sum(abs(A(i,[i+2:n])));
       if i==1;
           g(i)=A(i,2)/A(i,1);
           xs(i)=B(i)/A(i,1);
       else
           g(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*g(i-1));
           xs(i)=(B(i)-A(i,i-1)*xs(i-1))/(A(i,i)-A(i,i-1)*g(i-1));
       end
    elseif i>(m-2)
        otd=sum(abs(A(i,[1:i-2])));
        if i~=m;
        g(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*g(i-1));
        xs(i)=(B(i)-A(i,i-1)*xs(i-1))/(A(i,i)-A(i,i-1)*g(i-1));
        else
            g(i)=0;
            xs(i)=(B(i)-A(i,i-1)*xs(i-1))/(A(i,i)-A(i,i-1)*g(i-1));
        end
    else
        otd=sum(abs(A(i,[1:i-2,i+2:n])));
        g(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*g(i-1));
        xs(i)=(B(i)-A(i,i-1)*xs(i-1))/(A(i,i)-A(i,i-1)*g(i-1));
    end
    
    if otd~=0;
        display('A matrix is not Tridiagonal');
        x=[];        xs=[];        g=[];
        return;
    end
end

%% solve for x

x(m)=xs(m);
for i=(m-1):-1:1;
    x(i)=xs(i)-g(i)*x(i+1);
end

end