
% mex cpp (Only for the first time)
mex CMC_V0.cpp


% generate random tensor
size_tensor = [100,200,300];
tensor = double(rand(size_tensor)>0.7);

% CMC
rwu = CMC(tensor);

% extract exp(r), exp(w), exp(u)
r = result{1};
w = result{2};
u = result{3};

% probability tensor conditional on r, w, & u 
w = w';
u = reshape(u,[1,1,numel(u)]);
A = r.*w.*u;
P = A./(A+1);


