function regOpts = loadRegPresets(regFlag)

% you can modify the parameters after presetting to improve regularization
% performance e.g. "opts.R.p = 1.5;"

% the power of the gradient in norm
regOpts.p = 2;
regOpts.alphaTV = 0;
%opts.R.epsilon = 2*10^-8; % W12  
regOpts.alpha = 0;
regOpts.betta = 0;
regOpts.epsilon = 2*10^-5;
regOpts.gFlag = 0;

% 1 - use norms of the gradient of the model (Sobolev space norms), 
% 0 - use norms of the model itself
%opts.R.gFlag = 0;

switch regFlag
    case -2 % Tikhonov with derivative regularizaton
        regOpts.alpha = 100;
        regOpts.gFlag = 1;
        regOpts.alphaTV = 0;
        regOpts.betta = 0;
        regOpts.epsilon = 2*10^-5;
    case -1 % minimum support
        regOpts.betta = 5*10^-6; % Min Supp
        regOpts.epsilon = 10^-6; % Min Supp
        regOpts.gFlag = 0;
        regOpts.alphaTV = 0;
    case 1 % MGS
        regOpts.gFlag = 1;
        regOpts.epsilon = 1*10^-7; % MGS
        regOpts.betta = 5*10^-5; % MGS, TV
        regOpts.alphaTV = 0;
    case 2 % TV
        regOpts.betta = 0; % TV
        regOpts.gFlag = 1;
        regOpts.alphaTV = 0.01;
        regOpts.epsilon = 10^-8;
        regOpts.p = 1; % W12
    case 3 % W_p^1 p = 1.2
        regOpts.betta = 0; % TV
        regOpts.gFlag = 1;
        regOpts.alphaTV = 0.01;
        regOpts.epsilon = 10^-7;
        regOpts.p = 1.1; % W12
end

end