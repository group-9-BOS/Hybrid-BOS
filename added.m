%This function adds the coarse displacement from cross-correlation and
%refined from optical flow algorithm

function[u_final,v_final]=added(u_of,v_of,u_cc,v_cc)

addedu= u_of+u_cc;
addedv=v_of+v_cc;

u_final=addedu;
v_final=addedv;
end


























% load("ux.mat")
% dux=ux;
% load("uy.mat")
% duy=uy;
% load("quivUf.mat")
% ux=quivUf;
% load("quivVf.mat")
% uy=quivVf;
% addedu= dux+ux;
% addedv=duy+uy;
% figure;
% imagesc(addedv)
% 
% colorbar
% figure ;
% imagesc(addedu)
% 
% colorbar
