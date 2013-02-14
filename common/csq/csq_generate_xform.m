function [psi invpsi] = csq_generate_xform(xform,params)
% [psi inpsi] = csq_generate_xform(xform,param)
% Given a transform name and a set of accompanying parameters, this
% function will return two function handles which represent the forward and
% inverse transform.


switch xform
%     case 'dwt2d' 
%        [psi invpsi] = xform_dwt2d(params);
    otherwise
        return_str = sprintf('Transform "%s" is unsupported.',xform);
        error('csq_generate_xform:UnsupportedTransform',return_str);
end




%--------------------------------------------------

% function [psi invpsi] = xform_dwt2d(params)


