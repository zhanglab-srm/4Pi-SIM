function img_read=myimreadstack_TIRF(filename,nk,zhens,testreadx,testready)
    warning off;
	t = Tiff(filename, 'r');
    t.setDirectory(nk);  % start frame
    img_read = zeros(testreadx,testready,zhens);
    for k = 1:zhens-1
		img_read(:,:,k) = t.read();
        t.nextDirectory();
    end
    img_read(:,:,k+1) = t.read();
	t.close();
    warning on;
end