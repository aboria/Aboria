N = 100^2 %size of matrix
A = ones(N,N);
v = ones(N,1);
ncores = 2;
for i=1:ncores
    parpool(i)
    time = timeit(@() multiply(A,v));
    result(i,1) = i
    result(i,2) = time
    delete(gcp)
end

csvwrite('matlab_matrix_multiply_scaling.csv',result);

function multiply(A,v)
v = v + A*v;
end