function [A, rows,cols] = showme(filename, sqr)

% shows the matrix from an .mm file using imshow


[A, rows, cols] = mmread(filename);
%if sqr == 1, reshape to a square
if sqr==1
    dim = (rows*cols)^0.5;
    A = reshape(A, dim, dim);
    
end

rows
cols
A=full(A);
imshow(A);