function img_read=imreadstack_TIRF(filename,varargin)
    warning off;
	t = Tiff(filename, 'r');
    if numel(varargin)==0
        info = imfinfo(filename);
        t.setDirectory(1);
        num_images = numel(info);
        testreadx=info(1).Height;
        testready=info(1).Width;
    elseif numel(varargin)==1
        info = imfinfo(filename);
        t.setDirectory(varargin{1});
        num_images = numel(info)-varargin{1}+1;
        testreadx=info(1).Height;
        testready=info(1).Width;
    elseif numel(varargin)==2
        info = imfinfo(filename);
        t.setDirectory(varargin{1});
        num_images = varargin{2}-varargin{1}+1;
        testreadx=info(1).Height;
        testready=info(1).Width;
    elseif numel(varargin)==4  
        t.setDirectory(varargin{1});
        num_images = varargin{2}-varargin{1}+1;
        testreadx= varargin{3};
        testready= varargin{4};
    end
    img_read=zeros(testreadx,testready,num_images);
    for k = 1:num_images-1
		img_read(:,:,k)=t.read();
        t.nextDirectory();
    end
    img_read(:,:,num_images)=t.read();
	t.close();
    warning on;
end