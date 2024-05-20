function [n2,V,SA,x,y,z,chi2]=Chan_edge_detection_algorithm_Tong_modified(folder)
%% LOCATE POINTS DEFINING VACUOLE SURFACE FROM CONFOCAL Z-STACK

%intialize variables
%
% clear;

%figure
%hold on

xyvox=0.05;         %in microns, the xy pixel size after resampling
zvox=0.05;          %in microns, the z pixel size after resampling
xyvoxinit=0.065;    %in microns, the xy pixel size before resampling
zvoxinit=0.2;      %in microns, the z pixel size before resampling
%% Read in data images and manually-annotated center points

%folder=path with z-stack of images saved as individual TIFF files
% folder='/Users/tongshu/Desktop/Server download file/20201201-vacuole volume measurement/test3';

%read in data images from folder   
[N,M,S]=fileattrib(fullfile(folder,'*.tif'));
C=struct2cell(M);
names=C(1,1,:);
L=length(names);
for ntemp=1:L
    xtemp=im2double(imread(names{1,1,ntemp}));
    input(:,:,ntemp)=xtemp;
end
input=input-min(min(min(input))); %remove background
input=input/max(max(max(input))); %normalize pixel values from 0-1

%The folder should also contain a csv(text) file with a table of vacuole center
%point. We have typically done this by using ImageJ to select center
%points as ROIs using the multi-point selection tool. 

%%%Here the file is saved as csv file from imagej multi-point measurement.
%['filename = text file name in the quotes']
csvfile=dir([folder,'/*.csv']); 
csvfile=csvfile(~startsWith({csvfile.name},'.')); %remove csv file starting with dot

%read in the vacuole center coordinates from the csv file
init=importdata(fullfile(folder,csvfile.name)); %from csv importdata, 'init' is already a struct data.
% ij=table2struct(init);
ij = init;
    n2=size(ij.data,1); %Number of vacuoles counted from csv file
%     numCells=1;
%     indexVacs=1;
%     totalVacs=0;
    index_Y = strcmp(ij.colheaders,'Y');
    index_X = strcmp(ij.colheaders,'X');
    index_Slice = strcmp(ij.colheaders,'Slice');
    %Store vacolue pixel position in array
    x=round((1/xyvox)*ij.data(:,index_Y));
    y=round((1/xyvox)*ij.data(:,index_X));
    z=round((zvoxinit/zvox)*ij.data(:,index_Slice));            

clear N M S C L names xtemp; %clean up some variables

%% RESAMPLING INPUT

%This section resamples the z-stack between 2 cartesian coordinate systems:
%From: input voxel spacings defined by xyvoxinit and zvoxinit
%Into: output voxel spacings definedby xyvox and zvox

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %Original code pixel information position:
% % xyvox=0.05;         %in microns, the xy pixel size after resampling
% % zvox=0.05;          %in microns, the z pixel size after resampling
% % xyvoxinit=0.061;    %in microns, the xy pixel size before resampling
% % zvoxinit=0.15;      %in microns, the z pixel size before resampling
%I've found the actual z-spacing to generally be pretty different from the nominal
%spacing. More specifically, using the z-slice spacing from microscope settings
%tends to create reconstructed shapes that are noticeably elongated in z. 
%This seems to be an issue with index of refraction changes, and may be different
%system to system. I generally set this so that vacuoles aren't elongated
%or squashed.

%input=input-min(min(min(input)));
%inputres=myArrayResample3D([0.06 0.06 0.2],input,[0.06 0.06 0.06],'linear')
%
[rr,cc,dd] = size(input);
%r = r-1; c = c-1; d = d-1;

% smaller variable names
ops = [xyvoxinit xyvoxinit zvoxinit];
nps = [xyvox xyvox zvox];

% create the meshes for old and new image
%let meshes in x,y (ie where size is reduced) go from half-pix to size-halfpix
%let mesh in z (ie where size is increased) go from 1 to size
[Ox,Oy,Oz] = meshgrid(ops(2):ops(2):cc*ops(2),...
    ops(1):ops(1):rr*ops(1), ...
    ops(3):ops(3):dd*ops(3));
[Nx,Ny,Nz] = meshgrid(nps(2):nps(2):cc*ops(2),...
    nps(1):nps(1):rr*ops(1), ...
    nps(3):nps(3):dd*ops(3));%1:nps(3):d);

% create the new z-stack with name input with new cartesian spacing
input=interp3(Ox,Oy,Oz,input,Nx,Ny,Nz,'*cubic');
%}

contourxy(size(input,1),size(input,2),3)=0;

%%DRAWING RAYS AND FINDING INTENSITY MAX ALONG RAYS
%This section finds the boundary of each object specified in the input
%For each center point, _M_ rays are drawn emanating in random directions.
%Each ray is sampled at _N_ intervals of size _radian_ microns.

%N is the number of points along each ray drawn
%M is the number of angles (or equally spaced latitudinal angles, not usually used)
%(O is the number of equally spaced longitudinal angles, not usually used)

N=200;
M=500;
O=20;
M1=12;
radian=0.02;        %in microns, the sampling along each ray during intensity analysis

sizex=size(input,1);
sizey=size(input,2);
sizez=size(input,3);

%rs is the array of radii of the equally spaced points on the rays
%thetas is the array of latitudinal angles of each ray
%phis is the array of longitudinal angles of each ray.
%as given, the rand command makes a random set of rays that tend to be
%equally distributed around a sphere

rs=linspace(1,N+1,N+1);
thetas=2*pi*rand(M,1);%linspace(0,2*pi,M+1);%use linspace if you want equally spaced angles
phis=acos(2*rand(M,1)-1);%linspace(0,pi,O+1);%use phis if you want 3-D, don't if 2-D

%angles for plotting the ellipsoid fits in the end
thetasplot=linspace(0,2*pi,M1+1);
phisplot=linspace(0,pi,O+1);

%initialize variables
polx(N+1,M+1)=0;
poly(N+1,M+1)=0;
polz(N+1,M+1)=0;
polxfull(N+1,M+1)=0;
polyfull(N+1,M+1)=0;
polzfull(N+1,M+1)=0;
remx(N+1,M+1)=0;
remy(N+1,M+1)=0;
remz(N+1,M+1)=0;
C1(N+1,M+1,8)=0;
b=1;

%for each point on a ray (radian)
for i=1:N+1
    
    %for each ray (angle)
    for j=1:M
        
        %for each longitudinal angle
        %for k=1:O+1
            
            %calculate the rounded-down x,y,z-coordinates 'polx/y/z(i,j)' 
            %as well as the non-rounded coordinates 'polx/y/zfull(i,j)' 
            
            %polx(i,j,k)=int8(rs(i)*cos(thetas(j))*sin(phis(k)));
            %poly(i,j,k)=int8(rs(i)*sin(thetas(j))*sin(phis(k)));
            %polz(i,j,k)=int8(rs(i)*cos(phis(k)));
          
            polxfull(i,j)=(radian*rs(i)*cos(thetas(j))*sin(phis(j)))/xyvox;
            polyfull(i,j)=(radian*rs(i)*sin(thetas(j))*sin(phis(j)))/xyvox;
            polzfull(i,j)=(radian*rs(i)*cos(phis(j)))/zvox;
                     
            polx(i,j)=fix(polxfull(i,j));
            poly(i,j)=fix(polyfull(i,j));
            polz(i,j)=fix(polzfull(i,j));
            
            %remx/y/z(i,j) is the remainder from rounding-down, assembled
            %into the matrix C1 for sampling
            
            remx(i,j)=mod(polxfull(i,j)+1000,1);
            remy(i,j)=mod(polyfull(i,j)+1000,1);
            remz(i,j)=mod(polzfull(i,j)+1000,1);
            
            C1(i,j,8)=remx(i,j)*remy(i,j)*remz(i,j);
            C1(i,j,7)=remx(i,j)*remy(i,j)*(1-remz(i,j));
            C1(i,j,6)=remx(i,j)*(1-remy(i,j))*remz(i,j);
            C1(i,j,5)=remx(i,j)*(1-remy(i,j))*(1-remz(i,j));
            C1(i,j,4)=(1-remx(i,j))*remy(i,j)*remz(i,j);
            C1(i,j,3)=(1-remx(i,j))*remy(i,j)*(1-remz(i,j));
            C1(i,j,2)=(1-remx(i,j))*(1-remy(i,j))*remz(i,j);
            C1(i,j,1)=(1-remx(i,j))*(1-remy(i,j))*(1-remz(i,j));
    end
%    end
end


%for the number of center points
num=1; %Not sure about the meaning of num variable, since x,y,z are all colume arrays (n*1)
for m=1:n2
    %initialize variables
    f1{m}(N,M)=0; 
    diffinput=zeros(N,M);
    %for each point on the ray of index 2 and greater
    for i=2:N-1
        %for each ray (angle)
        for j=1:M
            %
            xp=polx(i+1,j)+x(m,num); %x coordinate of the index+1 point on the ray
            xc=polx(i,j)+x(m,num);   %x coordinate of the index   point on the ray
            xm=polx(i-1,j)+x(m,num); %x coordinate of the index-1 point on the ray
            yp=poly(i+1,j)+y(m,num); %y coordinate of the index+1 point on the ray
            yc=poly(i,j)+y(m,num);   %y coordinate of the index   point on the ray
            ym=poly(i-1,j)+y(m,num); %y coordinate of the index-1 point on the ray    
            zp=polz(i+1,j)+z(m,num); %z coordinate of the index+1 point on the ray
            zc=polz(i,j)+z(m,num);   %z coordinate of the index   point on the ray
            zm=polz(i-1,j)+z(m,num); %z coordinate of the index-1 point on the ray
            %
            
            %if the point with that radius and angle from the
            %center point (x(m,num), y(m,num)) is outside the original
            %image size, set the value to 0
            
            if xp+1>sizex || yp+1>sizey || zp+1>sizez || xp<1 || yp<1 || zp<1 || xm<1 || ym<1 || zm<1
                f1{m}(i,j)=0;

            %if the point is in the image, calculate the radial
            %graident at that point in the direction of the ray
            else
                %_inputc/p/m_ is the array of intensities in the 8-pixel
                %neighborhood around the pixel (xc/p/m,yc/p/m,zc/p/m)
                
                inputc=[input(xc,yc,zc) input(xc,yc,zc+1) input(xc,yc+1,zc) input(xc,yc+1,zc+1) input(xc+1,yc,zc) input(xc+1,yc,zc+1) input(xc+1,yc+1,zc) input(xc+1,yc+1,zc+1)];
                inputp=[input(xp,yp,zp) input(xp,yp,zp+1) input(xp,yp+1,zp) input(xp,yp+1,zp+1) input(xp+1,yp,zp) input(xp+1,yp,zp+1) input(xp+1,yp+1,zp) input(xp+1,yp+1,zp+1)];
                inputm=[input(xm,ym,zm) input(xm,ym,zm+1) input(xm,ym+1,zm) input(xm,ym+1,zm+1) input(xm+1,ym,zm) input(xm+1,ym,zm+1) input(xm+1,ym+1,zm) input(xm+1,ym+1,zm+1)];
                

                %}
                
                %_inputp/m/cinterp_ are the sum of intensities in the 8-pixel
                %neighborhood weighted by the remainders describing the
                %position of the points WRT the image lattice.
                
                inputpinterp=inputp*squeeze(C1(i+1,j,:));
                inputminterp=inputm*squeeze(C1(i-1,j,:));
                inputcinterp=inputc*squeeze(C1(i,j,:));
              
                %set _diffinput(i,j)_ to either the interpolated pixel
                %intensity (0th derivative, line 238) or the difference between the
                %i+1 and i-1 index pixel (~1st derivative, line 239)
                
                diffinput(i,j)=inputcinterp;
                %diffinput(i,j)=inputminterp-inputpinterp;

            end
        end
    end

    %set all values for the first/second/second-to-last/last radial point on all rays to 0
    f1{m}(1,:)=0; f1{m}(2,:)=0; f1{m}(N-1,:)=0; f1{m}(N,:)=0;
    
    %for other points along the ray
    for i=3:N-2
        for j=1:M
            %set deriv2 to be either the 0th derivative (interpolated
            %intensity) or ~2nd derivative (calculated by differences) at the point  
            deriv2=diffinput(i,j);
            
            %f1{m} is the collected values for each point _i_ on each ray
            %_j_. The value assigned is _deriv2_, possibly corrected by a
            %value to penalize for greater distances away from the origin
            %to make a preference for points closer vs. further away.
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   % This function weights the intensity values of pixels depending on
   % distance away from the origin. This is typically used to weight nearer
   % pixels more than distal points, in case of multiple membrane
   % crossings.
   
            f1{m}(i,j)=deriv2/(1+(i/N));%;%/i^0.5;%%/
            weight(i,1)=1;
        end
    end
    
    %C/Cmin are the max/min intensity (~2nd deriv) values along each ray,
    %I/Imin are the indices of those points along each ray.
    [C,I]=max(f1{m}(1:N-2,:),[],1);
    [Cmin,Imin]=min(f1{m}(:,:),[],1);
    s=std(I);       %standard deviation of the indices of the maxima
    avg=mean(I,2);  %mean of the indices of the maxima
    a=1;
    

    
    for j=1:M %for all rays

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %test to see if the index found for each ray is within a range determined by the
        %average (_avg_) and standard deviation (_s_) of intensities. This can reduce the
        %noise in the point cloud, but biases the values you get.
        if I(1,j)<=avg+2*s && I(1,j)>=avg-0.5*s %normally 2 0.5 %if I(1,j,k)<=avgk+sk && I(1,j,k)>=avgk-sk            
            %find the x,y,z coordinates associated with each index and ray
            x1=(polxfull(I(1,j),j)+x(m,1));
            y1=(polyfull(I(1,j),j)+y(m,1));
            z1=(polzfull(I(1,j),j)+z(m,1));
            
            if b>1 && x1==vertices(b-1,1) && y1==vertices(b-1,2) && z1==vertices(b-1,3)
            else
                vertices(b,:)=[x1*xyvox y1*xyvox zvox*z1 m];
                b=b+1;
                
                %assemble pixel coordinates of identified boundary points in _contourxy_    
                if round(x1)>0 && round(y1)>0 && round(z1)>0
                    contourxy(round(x1),round(y1),3,round(z1))=1;
                end

            end
            
            %assemble pixel positions in real space for fitting with
            %ellipsoids
            if a>1 && x1==vert{m}(a-1,1) && y1==vert{m}(a-1,2) && z1==vert{m}(a-1,3)
            else
                vert{m}(a,:)=[xyvox*x1 xyvox*y1 zvox*z1 thetas(j) phis(j)];
                a=a+1;
            end
        end
    end
end

%% FITTING THE POINT CLOUD WITH A SURFACE

clear p t tnorm V SA chi2;
W=500;

thetasout=2*pi*rand(W,1);
phisout=acos(2*rand(W,1)-1);

figure
for zz=1:n2
    clear X rho i j k ijk B cijk rhoout Xout x0 y0 z0 Xcent rho Bout phis thetas imax jmax kmax;
    clear cylX cylY cylZ p;
    
    %ellipsoid fitting
    %{
    thetasout=2*pi*rand(W,1);
    phisout=acos(2*rand(W,1)-1);
    
    thetas=thetasout;
    phis=phisout;
    
    X(:,1)=1*sin(phis).*cos(thetas)+rand(1);
    X(:,2)=1*sin(phis).*sin(thetas)+rand(1);
    X(:,3)=1*cos(phis)+rand(1);
    %}
    
    %data from edges and convert to center-of-mass coordinates
    X=vert{zz}(:,1:3);
    
    x0=mean(X(:,1));
    y0=mean(X(:,2));
    z0=mean(X(:,3));
    
    Xcent(:,1)=X(:,1)-x0;
    Xcent(:,2)=X(:,2)-y0;
    Xcent(:,3)=X(:,3)-z0;
        
    rho=((Xcent(:,1).^2+Xcent(:,2).^2+Xcent(:,3).^2).^0.5);
    %
    %for cylinder
    thetas=atan2(Xcent(:,2),Xcent(:,1));
    phis=acos(Xcent(:,3)./rho);
    %}
  
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%i/j/kmax determine how many basis functions are used to fit the point
%cloud. At least one of these should be less than or equal to 2 (usually kmax), the other
%two can be arbitrarily high. 3,3,2 are the values suggested in the
%reference.
    imax(zz)=3;
    jmax(zz)=3;
    kmax(zz)=2;

    %chi2 is a metric for goodness of the shape fit to the point cloud
    chi2(zz)=0;

    
  
    for ii=1:imax(zz)      
        i=ii-1;
        for jj=1:jmax(zz)
            j=(jj-1);
            for kk=1:kmax(zz)
                k=(kk-1);
                ijk=ii+j*imax(zz)+k*imax(zz)*jmax(zz);
                for nn=1:size(thetasout,1)
                    Bout(nn,ijk)=(sin(phisout(nn)))^(i+j)*(cos(phisout(nn)))^k*(sin(thetasout(nn)))^j*(cos(thetasout(nn)))^i;
                end
                for n=1:size(rho,1)
                    B(n,ijk)=(sin(phis(n)))^(i+j)*(cos(phis(n)))^k*(sin(thetas(n)))^j*(cos(thetas(n)))^i;
                end
            end
        end
    end
    cijk=B\rho;
    rhoout=Bout*cijk;
    rhochi=B*cijk;
    Xout=[rhoout.*cos(thetasout).*sin(phisout)+x0 rhoout.*sin(thetasout).*sin(phisout)+y0 rhoout.*cos(phisout)+z0];
   

    chi2(zz)=sum((rho-rhochi).^2);

%create triangle mesh to calculated surface fit using MyRobustCrust.
%MyRobustCrust calculates volume and surface area as well.
    p{zz}=[Xout(:,1) Xout(:,2) Xout(:,3)];
    [t{zz},tnorm{zz}]=MyRobustCrust(p{zz});
    V(zz)=SurfaceVolume(p{zz},t{zz},tnorm{zz});

    SA(zz)=0;
    for g=1:size(t{zz},1)
        aa=((p{zz}(t{zz}(g,1),1)-p{zz}(t{zz}(g,2),1))^2+(p{zz}(t{zz}(g,1),2)-p{zz}(t{zz}(g,2),2))^2+(p{zz}(t{zz}(g,1),3)-p{zz}(t{zz}(g,2),3))^2)^0.5;
        bb=((p{zz}(t{zz}(g,3),1)-p{zz}(t{zz}(g,2),1))^2+(p{zz}(t{zz}(g,3),2)-p{zz}(t{zz}(g,2),2))^2+(p{zz}(t{zz}(g,3),3)-p{zz}(t{zz}(g,2),3))^2)^0.5;
        cc=((p{zz}(t{zz}(g,3),1)-p{zz}(t{zz}(g,1),1))^2+(p{zz}(t{zz}(g,3),2)-p{zz}(t{zz}(g,1),2))^2+(p{zz}(t{zz}(g,3),3)-p{zz}(t{zz}(g,1),3))^2)^0.5;
        s=(aa+bb+cc)/2;
        area=(s*(s-aa)*(s-bb)*(s-cc))^0.5;
        SA(zz)=SA(zz)+area;
    end
    
    sphericity(zz)=pi^(1/3)*((6*V(zz))^(2/3))/SA(zz);
       
%rotate colors of point cloud and fits for making MATLAB figure
    if mod(zz,5)==2
        color=[0 1 0];
    elseif mod(zz,5)==3
        color=[1 0 0];
    elseif mod(zz,5)==4
        color=[0 1 1];
    elseif mod(zz,5)==0
        color=[1 1 0];
    elseif mod(zz,5)==1
        color=[1 0 1];
    end
    
  
    hold on   
    scatter3(X(:,1),X(:,2),X(:,3),16,'filled','MarkerFaceColor',color)

    
    
    trisurf(t{zz},p{zz}(:,1),p{zz}(:,2),p{zz}(:,3),'facecolor','black','edgecolor','black','FaceAlpha',0.2)
    axis equal
   

end


%write out volume, surface area, coordinates of "center", and chi2 value to
%text file
%{
fid=fopen('120209 vac8Dlog calcofluor z=140nm.txt','a');
for zz=1:n2
    fprintf(fid,'%s %1.0f %5.2f %5.2f %5.0f %5.0f %5.0f %10.5f\n',strcat(sample,'_',int2str(folderyes(folindex))),zz,V(zz), SA(zz), x(zz),y(zz),z(zz),chi2(zz));
end
%}
