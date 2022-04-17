
Rnn=0.1*eye(8);
Rgg=eye(7);
theta=[40;60;80;100;120;130;150];
vn=cos(theta*pi/180);
%----------------------%
%Computing drive matrix%
%----------------------%
for rows=1:8
    for columns=1:size(theta,1)
        A(rows,columns)=exp(1j*pi*(rows-1)*vn(columns));
    end

end

Rxx=A*Rgg*ctranspose(A)+Rnn;
[V,D]=eig(Rxx);
min=inf;
%------------------------------%
%Finding the minimum eigenvalue%
%------------------------------%
for i=1:8
    for j=1:8
    if j==i
      if min>=D(j,i)
          min=D(j,i);
          rows=j;
          columns=i;
      end
    end
    end
end
u=V(:,columns);%eigenvector of the min eigenvalue

%---------------------------------------------%
%Computing and plotting spatial power spectrum%
%---------------------------------------------%
angleobs=35;    
x=1;

max=-inf;
while angleobs<=151
    for rows=1:8
    ad(rows,1)=exp(1j*pi*cos(angleobs*pi/180)*(rows-1));
    end
    P(x,1)=1/((ctranspose(ad)*u)*ctranspose(u)*ad);
    angleobs=angleobs+0.01;
    
    if real(P(x,1))>=max
        max=real(P(x,1));
    end
    
     x=x+1;
end
P(:,1)=10*log10(P(:,1)/max);
plot(35:0.01:150.99,(real(P(:,1))));
