[N,M,S]=fileattrib(fullfile(subFolName,'*.tif'));



%M=dir(fullfile(folder,'*.tif'));
C=struct2cell(M);
names=C(1,1,:);

L=length(names);

for n=1:L
    xread=im2double(imread(names{1,1,n}));
    input(:,:,n)=xread;

end

input=input-min(min(min(input)));
input=input/max(max(max(input)));

%clear N M S C x L n names;
[Nfull,Mfull,Sfull]=fileattrib(fullfile(subFolName,'*.tif'));
Cfull=struct2cell(Mfull);
namesFull=Cfull(1,1,:);

Lfull=length(namesFull);
for n=Lfull+1:Lfull
    xread=im2double(imread(namesFull{1,1,n}));
    input(:,:,n)=xread;

end
         