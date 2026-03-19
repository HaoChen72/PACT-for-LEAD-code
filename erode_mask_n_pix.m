function mask = erode_mask_n_pix(mask, n)

    for i = 1:n 
        mask = imerode(mask, [0 1 0; 1 1 1; 0 0 0]);
    end

end