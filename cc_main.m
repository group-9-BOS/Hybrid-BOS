%This function implements cross-correlation algorithm for coarse
%displacement field

function[ux_cc,uy_cc]= cc_main(I1, I2)
wS     = 16;                    % Window size
sS     = 32;                    % Search size
intfac = 5;                     %this is the interpolation factor
th = 5; 
[ry cx] = size(I1);
rect = [1 1 size(I2,2) size(I2,1)];                                     % here we have set the region of interest to the entire region, here 1, 1  sets the coordinate of the top left corner of the image

nR_I = size(I1,1);                                                     % Number of rows in initial image
nC_I = size(I1,2);                                                     % Number of cols in initial image


wS2 = floor(wS/2);                                                    % Half of window size
sS2 = floor(sS/2);                                                    % Half of search size

%overlap of half of the window size is define basically 16:16:(500-8)=492

wIndR = (wS:wS:(nR_I-wS2))';                                  % Window center locations for each row 
wIndC = (wS:wS:(nC_I-wS2))';                                  % Window center locations for each column
numR  = length(wIndR);                                                      % Number of rows length returns the size of array
numC  = length(wIndC);                                                      % Number of columns returns the length of array


cPeak   = zeros(numR,numC);                                               % Initialize colPeak first initialized as zero
rPeak   = zeros(numR,numC);                                               % Initialize rowPeak the same
colOffset = zeros(numR,numC);                                               % Initialize colOffset the same
rowOffset = zeros(numR,numC);                                               % Initialize rowOffset the same


for i = 1:1:numR                                                            % Loop over all image rows
    for j = 1:1:numC                                                        % Loop over all image columns
        
        % Get window centers
        rCenter = wIndR(i);                                               % Window row center
        cCenter = wIndC(j);                                               % Window column center


        cropI1 = [cCenter-wS2, rCenter-wS2, wS, wS];        % Define cropping rectangle for image 1
        I1_Sub = imcrop(I1,cropI1);                                         % Crop image 1 (template)
        I1_Sub_transposed =  I1_Sub';
        % Crop the template to SEARCH size
        % - [col start, row start, col nums, row nums]
        cropI2 = [cCenter-sS2, rCenter-sS2, sS, sS];        % Define cropping rectangle for image 2
        I2_Sub = imcrop(I2,cropI2);                                         % Crop image 2 (comparison)
        I2_Sub_transposed =  I2_Sub';

 
        %here the cropI1 is the smaller image in size and this is then
        %moved along the cropI2 to get correlation. Crop I1 is 17x17 while
        % CropI2 is 33x33 because we started at (1,1)

% if part is not very necessary only kept as a caution 

 if (all(I1_Sub == I1_Sub(1,1)))                                     % If "window" array values are all the same
            rowP = 0;                                                       % No row offset
            colP = 0;                                                       % No column offset
            dx   = 0;
            dy   = 0;
        else
            % Compute normalized cross-correlation
            c = normxcorr2(I1_Sub,I2_Sub);                                  % Normalized cross correlation between the two cropped images 'c' gives us an array of values obtained by the correlation and the maximum is our actual coordinate
%             c = fftshift(fftshift(real(ifft2(conj(fft2(I1_sub).*fft2(I2_sub))), 1), 2); 
            [cR,cC] = size(c);                                              % cR no. of row in correlation map and cC no. of colmn in correlation map
            
            % Find the peak indices of the cross-correlation map
            [rowP,colP] = find(c == max(c(:)));                             % the value of location of Maximum of cross-correlation map is stored in the rowP and colP
        	rowP = rowP(1);                                                 % row where maximum occured 
            colP = colP(1);                                                 % column where the maximum occured
            

            %this is done simply so that array doesn't become out of bound
            % basically saying if the peak ours at first index take next
            % index and if the peak occurs at the last index take the
            % penultimate index
            
            
                if (rowP == 1)  rowP = rowP + 1; end      
                if (rowP == cR) rowP = rowP - 1; end
                if (colP == 1)  colP = colP + 1; end
                if (colP == cC) colP = colP - 1; end
                
% gaussian sub pixel interpolation
 
                dx   = (log(c(rowP-1,colP)) - log(c(rowP+1,colP)))/(2*log(c(rowP-1,colP)) - 4*log(c(rowP,colP)) + 2*log(c(rowP+1,colP)));
                dy   = (log(c(rowP,colP-1)) - log(c(rowP,colP+1)))/(2*log(c(rowP,colP-1)) - 4*log(c(rowP,colP)) + 2*log(c(rowP,colP+1)));

 % Set the col and row peak values from max of cross-correlation
        cPeak(i,j) = colP + dy;                                           % Column peak location (X)
        rPeak(i,j) = rowP + dx;                                           % Row peak location (Y)
        
        % Find the pixel offsets for X and Y directions- done to convert the
        % location from local to global.
        
        colOffset(i,j) = cPeak(i,j) - wS2 - sS2 - 1;                % Actual column pixel shift
        rowOffset(i,j) = rPeak(i,j) - wS2 - sS2 - 1;                % Actual row pixel shift


end
    end
end


% Meshgrid
[RR,CC] = meshgrid(wIndR,wIndC);                                            
RR      = RR';                                                              
CC      = CC';                                                              
quivX   = CC;                                                               
quivY   = RR;                                                               
quivU   = colOffset;                                                        
quivV   = rowOffset;                                                        
quivVel = sqrt(quivU.^2 + quivV.^2);                                        

quivU(:,1) = quivU(:,1) + 1;                                                % adds 1 to each element of first column
quivV(1,:) = quivV(1,:) + 1;                                                % adds 1 to each element to first row of matrix both done because of coordinate starts from 1,1 


XX     = quivX - min(min(quivX));
scaleX = rect(3)/max(max(XX));
XX     = XX*scaleX + rect(1);
YY     = quivY - min(min(quivY));
scaleY = rect(4)/max(max(YY));
YY     = YY*scaleY + rect(2);
quivU   = real(quivU);
quivV   = real(quivV);

    XXPlot = real(XX);
    YYPlot = real(YY);
    ZZPlot = real(quivVel);
    
    quivU(abs(quivU) > th) = nan;
quivV(abs(quivV) > th) = nan;
quivVel = sqrt(quivU.^2 + quivV.^2);

numXXInterp = intfac*size(XX,2);
    numYYInterp = intfac*size(YY,1);
    XSmooth = linspace(min(XX(1,:)),max(XX(1,:)),numXXInterp);
    YSmooth = linspace(min(YY(:,1)),max(YY(:,1)),numYYInterp);


%     
    
[XXPlot,YYPlot] = meshgrid(XSmooth,YSmooth);
    ZZPlot = interp2(XX,YY,quivVel,XXPlot,YYPlot, 'cubic');   
   
%% interpolation for converting the final values to the image size
    
quivX1 = quivX; 
[xX, yX] = meshgrid(1:size(quivX1, 2), 1:size(quivX1, 1)); %orginal elements

[xq, yq] = meshgrid(linspace(1, size(quivX1, 2), cx), linspace(1, size(quivX1, 1), ry)); %final coordinates
quivXf = interp2(xX, yX, quivX1, xq, yq, 'cubic'); %interpolation scheme

quivY1= quivY;  
[xY, yY] = meshgrid(1:size(quivY1, 2), 1:size(quivY1, 1));
[xqY, yqY] = meshgrid(linspace(1, size(quivY1, 2), cx), linspace(1, size(quivY1, 1), ry));
quivYf = interp2(xY, yY, quivY1, xqY, yqY, 'cubic');

quivU1 = quivU;  

[xU, yU] = meshgrid(1:size(quivU1, 2), 1:size(quivU1, 1));
[xqU, yqU] = meshgrid(linspace(1, size(quivU1, 2), cx), linspace(1, size(quivU1, 1), ry));
quivUf = interp2(xU, yU, quivU1, xqU, yqU, 'cubic');

quivV1 = quivV;  

[xV, yV] = meshgrid(1:size(quivV1, 2), 1:size(quivV1, 1));
[xqV, yqV] = meshgrid(linspace(1, size(quivV1, 2), cx), linspace(1, size(quivV1, 1), ry));
quivVf = interp2(xV, yV, quivV1, xqV, yqV, 'cubic');


quivVel1 = quivVel;  
[xVel, yVel] = meshgrid(1:size(quivVel1, 2), 1:size(quivVel1, 1));
[xqVel, yqVel] = meshgrid(linspace(1, size(quivVel1, 2), cx), linspace(1, size(quivVel1, 1), ry));
quivVelf = interp2(xVel, yVel, quivVel1, xqVel, yqVel, 'cubic');

% this portion of resizing done by the liu shen why not use this because
% the matrix size is not preciesly 500x500
% 
% [n0,m0]=size(quivU);
% [n1,m1]=size(I1);

% scale=round((n1*m1/(n0*m0))^0.5);
% XXX = imresize(quivX,scale);
% YYY = imresize(quivY, scale);
% ux0=imresize(quivU,scale);
% uy0=imresize(quivV,scale);
% ux0 = quivUf;
% uy0 = quivVf;

ux_cc=quivUf;
uy_cc=quivVf;


%% shifting the original background image for further processing
  [Im1_shift_FD]=image_shifting_fd(quivUf,quivVf,I1,I2);
% % 
    [Im1_shift_bil]=image_shifting_bilinear(quivUf,quivVf,I1,I2); 
% % %
     outputFilenamebi = 'image_with_bilinear.tif';
         outputFilenameFD = 'image_with_FD.tif';
% % 
    imwrite(Im1_shift_bil, outputFilenamebi);
            imwrite(Im1_shift_FD, outputFilenameFD);

%% outputs

end


