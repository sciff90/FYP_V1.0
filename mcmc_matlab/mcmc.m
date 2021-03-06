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

accepted = 0;
flg=0;
count=1;
sigma = 1.0;
k=2;
while k<=N    
	xi = THETA(k-1,:) + sigma.*randn(1,2*(order+1));
    xi(1) = 1.0;
	ytest = filter(xi(order+2:2*(order+1)),xi(1:order+1),u);
    

	if max(abs(y-ytest))>elim
		THETA(k,:) = THETA(k-1,:);
	else
		THETA(k,:) = xi;
		accepted = accepted+1;		
			
    	end
    
    if(mod(count,1000)==0 && flg==0)
        if(accepted/1000>0.3)
            sigma = sigma*1.2;
        elseif(accepted/1000<0.25)
            sigma = sigma/1.2;
        else
            flg=1
            THETA(1,:)=THETA(k,:);
            k = 1
            sigma
        end
        count=0;    
        accepted=0;
    end
    count = count+1;
    k = k+1;
end

figure(1)
for ii=1:(order+1)*2;
	subplot(2,order+1,ii)
	hist(THETA(:,ii),100)
end
figure(2)
for ii=1:(order+1)*2;
	subplot(2,order+1,ii)
	plot(theta(:,ii))
	title('Parallel chain realisations') 
end
mean(theta)
theta_0'

theta_0'
