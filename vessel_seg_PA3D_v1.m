function VeMask = vessel_seg_PA3D_v1(A, varargin)
    
    if nargin > 1
        BH_opt = varargin{1};
        FG_opt = varargin{2};
        Bi_opt = varargin{3};
    else
        % Options for Bowler-hat transformation
        BH_opt.si = 30; % max diameter of vessel
        BH_opt.no = 20; % number of orientations
        BH_opt.gauss_sig = 0.6; % Sigma for denoising the BH transformed image

        % Options for Frangi filter
        FG_opt.BlackWhite=false;
        FG_opt.FrangiScaleRange=[1 7];
        FG_opt.FrangiBetaOne = 0.5;
        FG_opt.FrangiBetaTwo = 16;


        % Options for binarizing
        Bi_opt.nhood_size = [7 11];
        Bi_opt.sensitivity = 1;
    end
 
    %% Execution
    
    % Bolwer-hat
    bh_A = A;
    bh_A = (bh_A - min(bh_A(:))) / (max(bh_A(:)) - min(bh_A(:)));
    
    bhf_A = bowlerHat2D(bh_A, BH_opt.si, BH_opt.no);
    
    
    % Frangi filter
    fgf_A = FrangiFilter2D(bhf_A,FG_opt);

    fgf_A(fgf_A<5e-7) = 0;
    
    
    %% Create binarized image
    % 
    % BHMask = fgf_A > 5e-7;
    
    T = adaptthresh(fgf_A,Bi_opt.sensitivity, NeighborhoodSize = Bi_opt.nhood_size, Statistic= "mean");
    
    VeMask = imbinarize(fgf_A, T);


end