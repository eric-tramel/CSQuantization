function csq_reset_seed()

if csq_in_octave
    rand_seed = time; % Current UTC in seconds 
else
    rand_seed = java.lang.System.currentTimeMillis; % Current UTC in milliseconds
end

if csq_in_octave
    rand('state',mod(rand_seed,2^32));
    randn('state',mod(rand_seed,2^32));
else
    s = RandStream('mcg16807','Seed',mod(rand_seed,2^32));
    RandStream.setDefaultStream(s);
end