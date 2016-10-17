function [] = write_csr_matrix(A,filename);
% Write A in compressed sparse row format (zero-based) to filename
% Text file format: 
% - one number per line
% - First number is m = number of rows
% - Second number is n = number of columns
% - Next m+1 numbers are the row pointers
% - Next nnz numbers are the column indices
% - Next nnz numbers are the matrix values

fid = fopen(filename,'w');
m = size(A,1);
n = size(A,2);
fprintf(fid,'%d\n',m);
fprintf(fid,'%d\n',n);
counter = 0;
for i=1:m,
    fprintf(fid,'%d\n',counter);
    counter = counter + nnz(A(i,:));
end;
fprintf(fid,'%d\n',counter);
counter = 0;
for i=1:m,
    for k=find(A(i,:)),
        fprintf(fid,'%d\n',k-1);
    end;
end;
for i=1:m,
    for k=find(A(i,:)),
        fprintf(fid,'%.16e\n',full(A(i,k)));
    end;
end;
fclose(fid);

