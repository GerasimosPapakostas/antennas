Rnn=0.1*eye(8);
Rgg=eye(2);
mindif=180;
counter=0;
%------------------------------------------------------%
%For a range of angles between 0-180 and with step 0.01%  
%------------------------------------------------------%
for theta1=0:0.01:180
    
    theta2=180-theta1;
    counter=counter+1;
%--------------------%
%Compute drive matrix%
%--------------------%
for rows=1:8
      A(rows,1)=exp(1j*pi*(rows-1)*cos(theta1));
      A(rows,2)=exp(1j*pi*(rows-1)*cos(theta2));
 end
  Rxx=A*Rgg*ctranspose(A)+Rnn;
  %---------------------------%
  %Find the minimum eigenvalue%
  %---------------------------%
  [V,D]=eig(Rxx);
   min1=inf;%finding the min eigenvalue
     for i=1:8
       for j=1:8
        if j==i
          if min1>=D(j,i)
              min1=D(j,i);
              rows=j;
              columns=i;
          end
        end
       end
    end
    u=V(:,columns);
    
    %--------------------------%
    %Compute the power spectrum%
    %--------------------------%
    angleobs=0;    
    x=1;
    max1=-inf;
  while angleobs<=180
    for rows=1:8
    ad(rows,1)=exp(1j*pi*cos(angleobs*pi/180)*(rows-1));
    end
    P(x,1)=1/((ctranspose(ad)*u)*ctranspose(u)*ad);
    angleobs=angleobs+0.01;
    
    if real(P(x,1))>=max1
        max1=real(P(x,1));
    end
    
     x=x+1;
  end
  P(:,1)=10*log10(P(:,1)/max1);
  
  %-------------------------------------------%
  %Find the two topical maxes and their places%
  %-------------------------------------------%
  [B,I]= findpeaks(squeeze(real(P)));
  [B1,I1]=maxk(B,2);
  %-----------------------------------------------------%
  %Find the minimum difference where the two maxes occur%
  %-----------------------------------------------------%
  if abs(I(I1(1,1),1)/100-I(I1(2,1),1)/100-0.02)<=mindif
      
    max2pos=I(I1(1,1),1)/100-0.01;
    max3pos=I(I1(2,1))/100 -0.01;
    mindif=abs(max3pos-max2pos);%minimum difference
    
  end

end



