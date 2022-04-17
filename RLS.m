Q=100;
a=0.98;
A=ones(8,6);
g=zeros(6,1);
theta=[50;70;100;110;130;160];
%----------------------%
%Computing drive matrix%
%----------------------%
 for rows=2:8
     for columns=1:6
      A(rows,columns)=exp(1j*pi*(rows-1)*cos(theta(columns,1)*pi/180));%drive matrix
     end
     
 end
 %---------%
 %Algorithm%
 %---------%
for q=1:2:Q
    if q~=1
    p=cos(2*pi*q/10);
    h=(a^(-1)*(y)^(-1)*x);
    y=(q/(q-1)*(a^(-1)*y-a^(-1)*h*ctranspose(x)*y));
    
    wrls=wrls+h*(conj(p)-ctranspose(x)*wrls);
    else
      wrls=0;
      y=10^6*eye(8);
    end
    for rows=1:6
        if rows==1
          g(rows,1)=cos(2*pi*q/10);
        else
            g(rows,1)=normrnd(0,1);
        end
    end
    for rows=1:8
        n(rows,1)=normrnd(0,sqrt(0.01));
    end
    x=A*g+n;   
end


for columns=8:-1:2
  wrls(:,columns)=[];
end
display(wrls);



%-------------------------------%
%Computing the radiation diagram%
%-------------------------------%
x=0;
angleobs=0;%angle to compute the radiation plot
    while angleobs<=180
        for rows=1:8
          a(rows,1)=exp(1j*(rows-1)*pi*cos(angleobs*pi/180));
        end
        x=x+1;
        AF(x,1)=ctranspose(wrls)*a(:,1);
        angleobs=angleobs+0.1 ;
    end
    %--------------%
    %Normalizing AF%
    %--------------%
    maxAF=-inf;
    for k=1:size(AF,1)
        if abs(AF(k,1))>=maxAF
            maxAF=abs(AF(k,1));
        end
    end
    
    for b=1:size(AF,1)
    normalized_AF(b,1)=(abs(AF(b,1))/maxAF);
    end
    
   %------------------------------%
   %Plotting the radiation diagram%
   %------------------------------%
   plot(0:0.1:180,normalized_AF(:,1));
   
   
   
   
   %-------------------------%
   %Computing the differences%
   %-------------------------%
   [pks1,locs1]=findpeaks(-normalized_AF(:,1));%peaks near zero
   [pks,locs]=findpeaks(normalized_AF(:,1));%max peaks
    min=inf(size(pks,1),1);
    minimum=inf(6,1);
   
   
 for rows=1:6
    for j=1:size(locs,1)
        if rows==1
            if pks(j,1)==1
              positionmax=locs(j,1)/10-0.1;
            end
            
        else
          if abs(theta(rows,1)-locs1(j,1)/10-0.1)<=minimum(rows,1)
            minimum(rows,1)=abs(theta(rows,1)-locs1(j,1)/10-0.1);
            positionmin(rows,1)=locs1(j,1)/10-0.1;
          end
        end
    end   
 end
  for rows=1:6
    if rows==1
        differencetheta(rows,1)=abs(positionmax(rows,1)-theta(rows,1));
    else
        differencetheta(rows,1)=abs(positionmin(rows,1)-theta(rows,1));
    end
        
  end
    
       
