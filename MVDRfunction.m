function [min1,max1,mean1,SD1,min2,max2,mean2,SD2,min3,max3,mean3,SD3] = MVDRfunction(SNR,mindif)
%MVDRfunction 
% min1->minimum value of Ди0
% max1->maximum value of Ди0
% mean1->mean value of Ди0
% SD1->standard deviation of Ди0
% min2->minimum value of Ди1,Ди2
% max2->maximum value of Ди1,Ди2
% mean2->mean value of Ди1,Ди2
% SD2->standard deviation of Ди1,Ди2
% min3->minimum value of SINR
% max3->maximum value of SINR
% mean3->mean value of SINR
% SD3->standard deviation of SINR
% SNR->Snr to compute the above values
% mindif->minimum difference between the angles и0,и1,и2

a=30;
b=150;

%--------------------%
%Producing the angles%
%--------------------%
for rows=1:1000
  for columns=1:3
      triuds(rows,columns)=(b-a).*rand+a;
  end
  %--------------------------------------------------------%
  %Check if the limit for minimum distance conditions is met
  %--------------------------------------------------------%
  while abs(triuds(rows,1)-triuds(rows,2))<mindif || abs(triuds(rows,2)-triuds(rows,3))<mindif || abs(triuds(rows,1)-triuds(rows,3))<mindif
    for columns=1:3
      triuds(rows,columns)=(b-a).*rand+a;
    end
  end
  
end


        %------------------%
        %Computing wmv,SINR%
        %------------------%
for i=1:1000 %repeat for each threeset of angles
    a=ones(8,1);%drive vector
    Rgg=eye(3);%correlation matrix
    Pn=10^(-SNR/10);%power of noise
    Rnn=Pn*eye(8); %correlation matrix of noise signals
    Rgigi=eye(2);
    Ainter=ones(8,2);%matrix for the interpolation

      
          for rows=1:8
             for columns=1:3
                A(rows,columns)=exp(1j*pi*(rows-1)*cos(triuds(i,columns)*pi/180));%drive matrix
             end
             Ainter(rows,:)=A(rows,[2 3]);
          end
          Rxx=A*Rgg*(ctranspose(A))+Rnn;
          ad=A(:,1);%drive vector
          wmv(:,i)=inv(Rxx)*ad;
           
    Ruu=Ainter*Rgigi*ctranspose(Ainter)+Rnn;
    SINR(i,1)=ctranspose(wmv(:,i))*(Rxx-Ruu)*wmv(:,i)/(ctranspose(wmv(:,i))*Ruu*wmv(:,i));
    
    
    
end

   
    
    %-------------------------------%
    %Computing the radiation diagram%
    %-------------------------------%
    minimum1=inf(1000,2);
for i=1:1000
    x=0;
    angleobs=0;%angle to compute the radiation plot
    while angleobs<=180
        for rows=1:8
          a(rows,1)=exp(1j*(rows-1)*pi*cos(angleobs*pi/180));
        end
        x=x+1;
        AF(x,i)=ctranspose(wmv(:,i))*a(:,1);
        angleobs=angleobs+0.1;
    end
    
    maxAF=-inf;
    for k=1:size(AF,1)
      if abs(AF(k,i))>=maxAF
            maxAF=abs(AF(k,i));
      end
    end
    for b=1:size(AF,1)
      normalized_AF(b,i)=(abs(AF(b,i))/maxAF);   
    end
    
end
   

   %-------------------------%
   %Computing the differences%
   %-------------------------% 
   minimum=inf(1000,2);
   
  
for i=1:1000
   
   expectedAngle0=triuds(i,1);
   expectedAngle1=triuds(i,2);
   expectedAngle2=triuds(i,3);
   [pks,locs]=findpeaks(normalized_AF(:,i));%all peaks for the i threeset
   [pks1,locs1]=findpeaks(-normalized_AF(:,i));%peaks near 0 to find Ди1,Ди2 
   for j=1:size(locs,1)%scan all the angle of peaks
  
       if pks(j,1)==1
           
           positionmax(i,1)=locs(j,1)/10-0.1;
       end
         
         
        if abs(expectedAngle1-locs1(j,1)/10-0.1)<=minimum(i,1) %find the min differences
            minimum(i,1)=abs(expectedAngle1-locs1(j,1)/10-0.1);
            positionmin(i,1)=locs1(j,1)/10-0.1;
  
        end
        if abs(expectedAngle2-locs1(j,1)/10-0.1)<=minimum(i,2)
           minimum(i,2)=abs(expectedAngle2-locs(j,1)/10-0.1);
           positionmin(i,2)=locs1(j,1)/10-0.1;
        end
    end
   differencetheta0(i,1)=abs(positionmax(i,1)-expectedAngle0);
   differencetheta1(i,1)=abs(positionmin(i,1)-expectedAngle1);
   differencetheta2(i,1)=abs(positionmin(i,2)-expectedAngle2);
end
 

%-------------------------%
%Passing the values to txt%
%-------------------------%
fileID=fopen('AoAdev_SINR.txt','wt');
for i=1:1000
fprintf(fileID,'\t%f',triuds(i,1));
fprintf(fileID,'\t');
fprintf(fileID,'\t%f',triuds(i,2));
fprintf(fileID,'\t');
fprintf(fileID,'\t%f',triuds(i,3));
fprintf(fileID,'\t');
fprintf(fileID,'\t%f',differencetheta0(i,1));
fprintf(fileID,'\t');
fprintf(fileID,'\t%f',differencetheta1(i,1));
fprintf(fileID,'\t');
fprintf(fileID,'\t%f',differencetheta2(i,1));
fprintf(fileID,'\t');
fprintf(fileID,'\t%f',SINR(i,1));
fprintf(fileID,'\t');
fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('AoAdev_SINR.txt','r');
formatSpec = '%f';
FILE=fscanf(fileID,formatSpec);
fclose(fileID);
%----------------------------------%
%Computing the mins,maxs,means,Stds%
%----------------------------------%
sum=0;
sumcounter=0;
for i=4:7:7000
    sumcounter=sumcounter+1;
    sum(sumcounter,1)=FILE(i,1);%passing the values of ДИ0
end

%min,max,mean,std value of ДИ0
min1=round(min(sum),3);
max1=round(max(sum),3);
mean1=round(mean(sum),3);
SD1=round(std(sum),3);


 sum2counter=0;
for i=5:7:7000
      j=i+1;
      sum2counter=sum2counter+1;
      sum_2(sum2counter,1)=FILE(i,1);
      sum_2(sum2counter,2)=FILE(j,1);
end
sum2=zeros(2000,1);
%making a new matrix to compute mean,std
for i=1:2000
    if i<=1000
        sum2(i,1)=sum_2(i,1);
    else
        sum2(i,1)=sum_2(i-1000,2);
    end
end
%min,max,sd value of Ди1,Ди2

min2=round(min(sum2),3);
max2=round(max(sum2),3);
mean2=round(mean(sum2),3);
SD2=round(std(sum2),3);


sum3counter=0;
for i=7:7:7000
    sum3counter=sum3counter+1;
    sum3(sum3counter,1)=FILE(i,1);%passing the values of SINR
end


%min,max,sd value of SINR
min3=round(min(sum3),3);
max3=round(max(sum3),3);
mean3=round(mean(sum3),3);
SD3=round(std(sum3),3);



end

