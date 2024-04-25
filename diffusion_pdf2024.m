function y = diffusion_pdf2024(data,C,l1,l2,alpha_in)

data = sort(data);

integral_function = @(a, C, l1, l2, x) integral(@(t) ...
(C * l1 * l2 * ( (x - t).^a ) .* exp(C^2 ./ (2 * (t - x))) .* (exp(-l2 * t) - exp(-l1 * t))) / (sqrt(2 * pi) * (l1 - l2)), 0, x);

yy = arrayfun(@(x_val) integral_function(alpha_in, C, l1, l2, x_val),data);
yy(isnan(yy)) = 0;
yy(yy<0) = 0;

y = yy;