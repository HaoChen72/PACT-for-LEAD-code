function imda = bowlerHat2D(im,si,no)
  % im: 2D Gray-Scale Image
  % si: Expected Maximum Width of the Vessel
  % no: Number of Directions
  n_orien = 0:180/no:180-180/no;
  
  imol = zeros(size(im,1),size(im,2),length(si),length(n_orien));
  imod = zeros(size(im,1),size(im,2),length(si));
  for i=1:length(si)
    for j=1:length(n_orien)
      se = strel('line',si(i),n_orien(j));
      imol(:,:,i,j) = imopen(im,se);
    end
    se = strel('disk',round(si(i)/2));
    imod(:,:,i) = imopen(im,se);
  end
  
  imd = zeros(size(im,1),size(im,2),length(si));
  imr = zeros(size(im,1),size(im,2),length(si));
  imm = zeros(size(im,1),size(im,2),length(si));
  
  triv = imod==0;
  for i=1:length(si)
    imm(:,:,i) = max(squeeze(imol(:,:,i,:)),[],3);
    imd(:,:,i) = imm(:,:,i) - imod(:,:,i);
  end
  imr(triv) = 0;
  imda = max(imd,[],3);
  imda= double(imda); imda = (imda - min(imda(:))) / (max(imda(:)) - min(imda(:)));
end