i = 2.0;
ii = 0;
maxNumCompThreads(1)
while i<200
    ii = ii + 1;
    N = int32(i^2); %size of matrix
    N
    A = ones(N,N);
    v = ones(N,1);
    time = timeit(@() multiply(A,v))
    result_size(ii,1) = double(N);
    result_size(ii,2) = double(N*N)/time;
    i = i * 1.05; 
end

csvwrite('matlab_matrix_multiply_size.csv',result_size);

N = 100^2 %size of matrix
A = ones(N,N);
v = ones(N,1);
ncores = 32;
for i=1:ncores
    maxNumCompThreads(i)
    %parpool(i)
    time = timeit(@() multiply(A,v));
    result_scale(i,1) = i;
    result_scale(i,2) = N*N/time
    %delete(gcp)
end

csvwrite('matlab_matrix_multiply_scaling.csv',result_scale);

function multiply(A,v)
v = v + A*v;
end
