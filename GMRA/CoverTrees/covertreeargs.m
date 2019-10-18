function A=covertreeargs() 

%  if ne(strcmp(class(X),'double'),true)
%  'class(X) must be double'
%  A=[];
%  return
%end


%  A.X=X;

  A.theta=input('Input theta\n');
if ne(strcmp(class(A.theta),'double'),true)
'class of theta must be double'
A=[];
return
end

  A.numlevels=int32(input('Input numlevels\n'));
  A.minlevel=int32(input('Input minlevel\n'));
  A.NTHREADS=int32(input('Input NTHREADS\n'));
  if A.NTHREADS==0
    A.BLOCKSIZE=int32(0);
else
  A.BLOCKSIZE=int32(input('Input BLOCKSIZE\n'));
end

  return
