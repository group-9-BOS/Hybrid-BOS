function [Im1_shift_bil]=image_shifting_bilinear(ux,uy,I1,I2)


I1=double(I1);
I2=double(I2);
ux=double(ux); 
uy=double(uy);


[m0,n0]=size(ux);
[m1,n1]=size(I2); 


xB = ux;
yB = uy;

Im1_shift0=I2;%this us basically storing origial image
for i=1:m1
    for j=1:n1
        i_shift=i+round(yB(i,j));
        j_shift=j+round(xB(i,j));
        if (i_shift<=m1) && (i_shift >=1) && (j_shift<=n1) && (j_shift>=1)
            Im1_shift0(i_shift,j_shift)=I1(i,j);           
        else
            Im1_shift0(i,j)=I1(i,j);
        end
    end
end
Im3=Im1_shift0; %is I2 gone to I1
Im1_shift1=Im3;
ex=xB-round(xB); %error in x 
ey=yB-round(yB); %error in y

mask_size=10;
std=0.6*mask_size; %the standard deviation has to be found using hit and trial
H1=fspecial('gaussian',mask_size,std);
ex=imfilter(ex,H1);
ey=imfilter(ey,H1);


for i=1:(m1-1) %goes upto 1 less size of the image because we see that dx and dy =1 are being added in not arrayout of bound
    for j=1:(n1-1)
 N0 = (1-ex(i))*(1-ey(j));
 N1 = (1-ex(i))*ey(j);
 N2 = ex(i)*(1-ey(j));
 N3  =ex(i)*ey(j);
 
  Im1_shift1(i,j) = N0*Im1_shift0(i,j)+N1*Im1_shift0(i,j+1)+N2*Im1_shift0(i+1,j)+N3*Im1_shift0(i+1,j+1);

    end
end

Im1_shift_bil=uint8(Im1_shift1);


