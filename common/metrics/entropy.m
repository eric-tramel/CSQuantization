 function H = entropy(x,b)
% Calculate the entropy, empirically, of the values of x. This function
% assumes that the list of terms/values is limited to some dictionary of
% integer vales. Using any kind of fractional values will, more than
% likely, throw the entire calculation off.
% Sungkwang touched this file.

% Lets determine how many values we have in x
c = unique(x);
N = length(x);

% Now, lets determine the distribtuion of values in x
H = 0; 
for i=1:length(c)
    term = c(i);
    
    % Count up how many we have
    term_count = sum(x == term);
    
    % Now we get this distribution
    p_term = term_count ./ N;
    
    % Then, we get the amount of entropy added by this term
    H = H + p_term .* (log(p_term)./log(b));
end

H = -H;
