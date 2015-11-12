M = [1 2 3;1 2 3];
fid = fopen('output.csv','w');
fprintf(fid,'%s\r\n',['One,Two,Three']);
fclose(fid);
dlmwrite('output.csv', M,'-append','delimiter',',','roffset',0,'coffset',0);

N = csvread('output.csv',1,0)