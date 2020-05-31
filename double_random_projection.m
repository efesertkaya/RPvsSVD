% we will read images from subfolder images first get the names in order to
% read them

clear
clc
listing=dir('Images/lfwdataset');% read the imagenames

listing_cell=struct2cell(listing);% cell format

% ATTENTION DONOT RUN THE CODE DIRECTLY TO SEE RESULTS OF THE CODE OPEN THE
% FOLLOWING USING CTRL+O: \metric_results.mat 

imagelar=cell(length(listing_cell)-5,1);
for i=6:length(listing_cell)
imagelar{i-5,1}=listing_cell{1,i};
end

cd(listing_cell{2,1}) % to read the images change directory

times_dim=zeros(0,3);% first column svd comptime second column rp comp time third column no of dim 



it=1;


statstat=cell(51,2);

eps=0.1; % \epsilon used in JL lemma condition
kk= ceil((3/(96*eps^2))*log(3840));% how many singular values from matrix B_1


for iii=0:50

kk=kk+25*iii;
referenced_metrics=zeros(length(imagelar),3);% 1 psnr , 2 Mean squared error, 3 structural similarity  
non_referenced_metrics=zeros(length(imagelar),2);
referenced_metrics_rp=zeros(1000,3);% 1 psnr , 2 Mean squared error, 3 structural similarity , 4
non_referenced_metrics_rp=zeros(1000,2);


for j=1:length(imagelar)
A=imread([imagelar{j}]);
A=rgb2gray(A);
A=double(A); % convert data matrix to double to perform SVD
no_of_dim=kk;% number of first how much singular values to be selected
times_dim(it,3)=kk;

St_svd = tic;
[U,S,V]= svd(A); % SVD



Anew= U(:,1:no_of_dim) * S(1:no_of_dim,1:no_of_dim) * transpose(V(:,1:no_of_dim));

svd_time=toc(St_svd);
times_dim(it,1)=svd_time;

Anew=uint8(Anew);
imwrite(Anew,[listing_cell{2,1} '\001SVDVersions\' num2str(j) '_grSc_SVD_n_dim_' num2str(kk) '.jpg'],'jpg')

referenced_metrics(j,1)=psnr(Anew,uint8(A));
referenced_metrics(j,2)=immse(Anew,uint8(A));
referenced_metrics(j,3)=ssim(Anew,uint8(A));


non_referenced_metrics(j,1)=brisque(Anew);
non_referenced_metrics(j,2)=niqe(Anew);


start_RP=tic;

P11=normrnd(0,1,[size(A,1),kk]);

B11= (1/sqrt(kk))* P11'* A;

[U11,S11,V11]=svd(B11); 

hatV11=V11(:,1:no_of_dim);

Anewrp=A*(hatV11)*hatV11';

rp_time=toc(start_RP);
times_dim(it,2)=rp_time;

Anewrp=uint8(Anewrp);
imwrite(Anewrp,[listing_cell{2,1} '\000RPVersions\' num2str(j) '_grSc_RP_n_dim_' num2str(kk) '.jpg'],'jpg')


referenced_metrics_rp(j,1)=psnr(Anewrp,uint8(A));
referenced_metrics_rp(j,2)=immse(Anewrp,uint8(A));
referenced_metrics_rp(j,3)=ssim(Anewrp,uint8(A));



non_referenced_metrics_rp(j,1)=brisque(Anewrp);
non_referenced_metrics_rp(j,2)=niqe(Anewrp);


A=uint8(A);
imwrite(A,[listing_cell{2,1} '\000grayScaleVersions\' num2str(j) '_grSc.jpg'],'jpg')

it=it+1;
end

statstat{iii+1,1}=[referenced_metrics non_referenced_metrics];
statstat{iii+1,2}=[referenced_metrics_rp non_referenced_metrics_rp];

end

kk=26;
avg_stat_for_dim=zeros(0,13);
for i=0:13
    i
    kk=kk+i*25;
    temp1=statstat{i+1,1};
    temp2=statstat{i+1,2};
    
    temp2=temp2(1:45,:);
    
    temp3=times_dim(times_dim(:,3)==kk,:);
    stat_row=[mean(temp1) mean(temp2) mean(temp3)];    
    
    avg_stat_for_dim=[avg_stat_for_dim; stat_row];
end

figure(1)
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,1))
xlabel('Number of Dimensions')
ylabel('Average PSNR Value')
hold on
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,6))
legend('SVD','RP')
% 1 psnr , 2 Mean squared error, 3 structural similarity  
figure(2)
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,2))
xlabel('Number of Dimensions')
ylabel('Average MSE Value')
hold on
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,7))
legend('SVD','RP')

figure(3)
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,3))
xlabel('Number of Dimensions')
ylabel('Average SSIM Value')
hold on
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,8))
legend('SVD','RP')

figure(4)
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,4))
xlabel('Number of Dimensions')
ylabel('Average Brisque Value')
hold on
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,9))
legend('SVD','RP')

figure(5)
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,5))
xlabel('Number of Dimensions')
ylabel('Average Niqe Value')
hold on
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,10))
legend('SVD','RP')

figure(6)
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,11))
xlabel('Number of Dimensions')
ylabel('Average Computation Time')
hold on
plot(avg_stat_for_dim(:,13),avg_stat_for_dim(:,12))
legend('SVD','RP')


A=[16     2     3    13; 5    11    10     8;9     7     6    12; 4    14    15     1];

Aref=[ 2     2     3    1;     5    1    6     8;     9     7     6    12;     4    14    15     1];

ssim(uint8(A),uint8(Aref))

% Singular Value Decomposition Image Compression 

% I=imread('kleo.jpeg'); % imagename
% 
% no_of_dim=30;% number of first how much singular values to be selected
% A1=I(:,:,1);
% 
% A1=double(A1);
% [U1, S1,V1]= svd(A1);
% 
% A1new= U1(:,1:no_of_dim) * S1(1:no_of_dim,1:no_of_dim) * transpose(V1(:,1:no_of_dim));
% 
% A2=I(:,:,2);
% 
% A2=double(A2);
% [U2, S2,V2]= svd(A2);
% 
% A2new= U2(:,1:no_of_dim) * S2(1:no_of_dim,1:no_of_dim) * transpose(V2(:,1:no_of_dim));
% 
% A3=I(:,:,3);
% 
% A3=double(A3);
% [U3, S3,V3]= svd(A3);
% 
% A3new= U3(:,1:no_of_dim) * S3(1:no_of_dim,1:no_of_dim) * transpose(V3(:,1:no_of_dim));
% 
% A3new=uint8(A3new);
% A2new=uint8(A2new);
% A1new=uint8(A1new); 
% Inew= zeros(size(I));
% Inew=uint8(Inew);
% 
% Inew(:,:,1)= A1new;
% Inew(:,:,2)=A2new;
% Inew(:,:,3)=A3new;
% 
% figure(1)
% imshow(Inew)
% xlabel('Singular Value Decomposition')
% 
% 
% % Double random projection Algorithm
% k1=100; % number of columns in the first projection matrix
% k2=500; % number of columns in the second projection matrix
% no_of_dim2=50;
% 
% P11=normrnd(0,1,[size(I,1),k1]);
% 
% B11= (1/sqrt(k1))* P11'* A1;
% 
% [U11,S11,V11]=svd(B11); 
% 
% hatV11=V11(:,1:no_of_dim2);
% 
% Ak1=A1*(hatV11)*hatV11';
% 
% P21=normrnd(0,1,[size(I,2),k2]);
% 
% B21= A1* P21;
% 
% A1newrp =B21* pinv(hatV11'*P21) * hatV11';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% P12=normrnd(0,1,[size(I,1),k1]);
% 
% B12= (1/sqrt(k1))* P11'* A1;
% 
% [U12,S12,V12]=svd(B12); 
% 
% hatV12=V12(:,1:no_of_dim2);
% 
% Ak2=A2*(hatV12)*hatV12';
% 
% P22=normrnd(0,1,[size(I,2),k2]);
% 
% B22= A2* P22;
% 
% A2newrp =B22* pinv(hatV12'*P22) * hatV12';
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% P13=normrnd(0,1,[size(I,1),k1]);
% 
% B13= (1/sqrt(k1))* P13'* A1;
% 
% [U13,S13,V13]=svd(B13); 
% 
% hatV13=V13(:,1:no_of_dim2);
% 
% Ak3=A3*(hatV13)*hatV13';
% 
% P23=normrnd(0,1,[size(I,2),k2]);
% 
% B23= A3* P23;
% 
% A3newrp =B23* pinv(hatV13'*P23) * hatV13';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A1newrp=uint8(A1newrp);
% A2newrp=uint8(A2newrp);
% A3newrp=uint8(A3newrp); 
% Inewrp= zeros(size(I));
% Inewrp=uint8(Inewrp);
% 
% Inewrp(:,:,1)= A1newrp;
% Inewrp(:,:,2)=A2newrp;
% Inewrp(:,:,3)=A3newrp;
% 
% figure(2)
% imshow(Inewrp)
% xlabel('Double Random Projection')






