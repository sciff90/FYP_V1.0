clear all;
N = 1e5;
num_samples = 50;
order = 1;
elim = 0.1;


u = ones(1,num_samples);
[b,a] = butter(order,0.2);
z = filter(b,a,u);
theta_0 = [a b];
e = elim*(2*rand(size(z))-1);
y = z+e;

THETA = zeros(N,2*(order+1));  THETA(1,:)=theta_0;

% Generate a proposal;

for k=2:N

  xi = THETA(k-1,:) + 0.001*randn(1,2*(order+1));
  ytest = filter(xi(order+2:2*(order+1)),xi(1:order+1),u);

  if max(abs(y-ytest))>elim
      THETA(k,:) = THETA(k-1,:);
  else
      THETA(k,:) = xi;
  end;


end;

figure(2)
hist(THETA(:,2),100);
figure(3)
hist(THETA(:,3),100);

