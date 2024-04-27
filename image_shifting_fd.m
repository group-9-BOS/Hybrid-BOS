function [Im1_shift_FD]=image_shifting_fd(ux,uy,I1,I2)



Im1=double(I1);
Im2=double(I2);
uxI=double(ux);
uyI=double(uy);

[m0,n0]=size(ux);
[m1,n1]=size(Im1);
% generate a shifted image from Im1 based on the velocity field that is
% rounded

%this loop is used for shifting image 2 towards the image one. Im1_shift0
%is after all I2 and Im1 is I1
Im1_shift0=Im2;
for i=1:m1
    for j=1:n1
        i_shift=i+round(uyI(i,j));
        j_shift=j+round(uxI(i,j));
        if (i_shift<=m1) && (i_shift >=1) && (j_shift<=n1) && (j_shift>=1)
            Im1_shift0(i_shift,j_shift)=Im1(i,j);           
        else
            Im1_shift0(i,j)=Im1(i,j);
        end
    end
end


Im3=Im1_shift0; %is I2
Im1_shift1=Im3;
duxI=uxI-round(uxI); %sub pixels values 
duyI=uyI-round(uyI); %sub pixel values 

mask_size=10;
std=0.6*mask_size;
H1=fspecial('gaussian',mask_size,std);
duxI=imfilter(duxI,H1);
duyI=imfilter(duyI,H1);


for i=1:(m1-1)
    for j=1:(n1-1)
          gradx(i,j)=(Im3(i,j+1)*duxI(i,j+1)-Im3(i,j)*duxI(i,j));
          grady(i,j)=(Im3(i+1,j)*duyI(i+1,j)-Im3(i,j)*duyI(i,j));   
          Im1_shift1(i,j)=Im3(i,j)-(gradx(i,j)+grady(i,j));
    end
end


Im1_shift_FD=uint8(Im1_shift1);





