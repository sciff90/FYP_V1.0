clear all;
N = 1e7;
num_samples = 300;
order = 1;
elim = 0.1;


u = randn(1,num_samples);
[b,a] = butter(order,0.2);
[b,a] = cheby1(order,0.9,0.2);
z = filter(b,a,u);
theta_0 = [a b]';
theta_0(1) = 1.0;
e = elim*(2*rand(size(z))-1);
y = z+e;

tic;
theta = mcmc(u,y,N,order,theta_0,elim);
toc;
figure(1)
for ii=1:(order+1)*2;
	subplot(2,order+1,ii)
	hist(theta(:,ii),100)
end
figure(2)
for ii=1:(order+1)*2;
	subplot(2,order+1,ii)
	plot(theta(1:12:N,ii))
	title('Single chain realisations') 
end
figure(3)
for ii=1:(order+1)*2;
	subplot(2,order+1,ii)
	plot(theta(:,ii))
	title('Parallel chain realisations') 
end
mean(theta)
theta_0'

