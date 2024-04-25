function y = diffusion_pdf_conv2024(data,mu1,mu2,mu3,mu4)

data = sort(data);

Fexp = @(x,l1,l2) (x>=0) .* ( (l1.*l2/(l1-l2)) .*(exp(-l2.*x)-exp(-l1.*x)) ); 
Fdiff = @(x,C,a) (x>=0) .* (C/sqrt(2*pi)) .* (x.^a) .* exp(-(C^2)./(2*x));

f1 = Fexp(data,mu2,mu3);
f2 = Fdiff(data,mu1,mu4);
y = conv(f1,f2,'same');

y(isnan(y)) = 0;
y(y<0) = 0;
