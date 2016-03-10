% Inverse Transform Sampling

function rand_num = gen_rand_pdf(glob_level, glob_pdf, rand_num_count)

valid_glob_pdf = glob_pdf(glob_pdf>0);
valid_glob_level = glob_level(glob_pdf>0);

d_level = abs(glob_level(2) - glob_level(1));
cdf = cumsum(valid_glob_pdf)*d_level;

u = rand([rand_num_count,1]);

rand_num = interp1(cdf, valid_glob_level, u);


