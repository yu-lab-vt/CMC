




%% 2D tensor
tensor_size_2D = [11980,8952];
x = (rand(tensor_size_2D)>0.7);

% Get margin totals
nDim = numel(size(x));
yA = [];
yA{1,1}=double(sum(x,2));
yA{2,1}=double(sum(x,1));

% store to .bin files
dirTmp='.\CMC_2D\margins\';
for iDim=1:nDim
    yTmp=yA{iDim,1};
    fileID = fopen([dirTmp,'m',num2str(iDim-1),'.bin'],'w');
    fwrite(fileID,yTmp,"double");
    fclose(fileID);
end






%% 3D tensor
tensor_size_3D = [100,200,300];
x = (rand(tensor_size_3D)>0.7);

% Get margin totals
nDim = numel(size(x));
yA = [];
tmp = double(sum(x,[2,3]));
yA{1,1}=tmp(:);
tmp = double(sum(x,[1,3]));
yA{2,1}=tmp(:);
tmp = double(sum(x,[1,2]));
yA{3,1}=tmp(:);


% store to .bin files
dirTmp='.\CMC_3D\margins\';
for iDim=1:nDim
    yTmp=yA{iDim,1};
    fileID = fopen([dirTmp,'m',num2str(iDim-1),'.bin'],'w');
    fwrite(fileID,yTmp,"double");
    fclose(fileID);
end


