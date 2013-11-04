clear all;
N = 20000000;
num_samples = 50;
order = 2;
elim = 0.01;


u = ones(1,num_samples);
[b,a] = butter(order,0.2);
z = filter(b,a,u);
theta_0 = [a b]';
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

