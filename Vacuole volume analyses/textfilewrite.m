fid=fopen('test.txt','at')%,'n','Macintosh')
for x=1:3
    fprintf(fid,'%s %1.0f\n','testing', x);
end
fclose(fid);
%fprintf(fid,'%s %1.0f %5.2f %5.2f %5.0f %5.0f %5.0f %10.5 %1.0f\n',strcat(sample,'_',int2str(folderyes(folindex))),zz,V(zz), SA(zz), x(zz),y(zz),z(zz),chi2(zz),ij(zz,1));
