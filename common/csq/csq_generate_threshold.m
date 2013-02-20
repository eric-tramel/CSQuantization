function thresh = csq_generate_threshold(name,params)
% thresh = csq_generate_threshold(name,params)
% Generate a function handle which implements the chosen
% thresholding strategy.

switch name
    case 'top'
        csq_required_parameters(params,'k');
        thresh = @(x) top_k(x,params.k);
    case 'bivariate-shrinkage'
        csq_required_parameters(params,'imsize','L','end_level','windowsize','lambda');
        thresh = @(x) csq_dwt_cell2vec(bivariate_shrinkage(csq_dwt_vec2cell(x,...
                                                                            params.imsize(1), ...
                                                                            params.imsize(2), ...
                                                                            params.L), ...
                                                           params.lambda,...
                                                           params.windowsize,...
                                                           params.end_level));
    otherwise
        return_str = sprintf('Threshold "%s" is unsupported.',name);
        error('csq_generate_threshold:UnsupportedTransform',return_str);
end