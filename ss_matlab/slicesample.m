clear all;
N = 1e5;
num_samples = 20;
order = 1;
elim = 0.1;


u = ones(1,num_samples);
[b,a] = butter(order,0.2);
z = filter(b,a,u);
theta_0 = [a b];
e = elim*(2*rand(size(z))-1);
y = z+e;

theta = zeros(N,2*(order+1));  theta(1,:)=theta_0;theta(2,:)=theta_0;

width = 0.01;
k=2;

theta_new = theta(k-1,:);
pstar = ptarget(y,u,elim,theta_new,order)
	
while(k<=N)

	Puprime = pstar*rand; 
    Puprime = pstar+log(rand);
	
	k
    ii=2;
	while( ii<=2*(order+1))
		
		theta_new = theta(k,:);
		theta_new_l = theta(k,:);	
		theta_new_r = theta(k,:);
		theta_prime = theta(k,:);
	
		bit = rand;
		theta_new_l(ii)	= theta_new(ii)-bit*width;
		theta_new_r(ii)	= theta_new(ii)+(1-bit)*width;
		
		% step out until horizontal span density
		while(ptarget(y,u,elim,theta_new_l,order)>Puprime)
			theta_new_l(ii) = theta_new_l(ii) -width;       
		end
		
		while(ptarget(y,u,elim,theta_new_r,order)>Puprime)
			theta_new_r(ii) = theta_new_r(ii) +width;
		end
		
		stepcount = 0;
		
		while(1)
			stepcount = stepcount+1;
			% Draw a candidte value uniformly distributed in interval [Mleft,Mright]
			theta_prime(ii) = rand()*(theta_new_r(ii) - theta_new_l(ii)) + theta_new_l(ii);
                       
			pstar = ptarget(y,u,elim,theta_prime,order); 
            
			
			if(pstar>Puprime)                
				break;
			else
				if(theta_prime(ii)>theta_new(ii))
					theta_new_r(ii) = theta_prime(ii);
				elseif(theta_prime(ii)<theta_new(ii))
					theta_new_l(ii) = theta_prime(ii);
                else           
                    fprintf('ERROR\n');
                    break;
				end
			end					
        end
        theta(k,ii) = theta_prime(ii);
        ii=ii+1;		
    end
    theta(k,1) = 1.0;
    k=k+1;
    theta(k,:) = theta(k-1,:);

end

figure(1)
for ii=1:(order+1)*2;
	subplot(2,order+1,ii)
	hist(theta(:,ii),100)
end;
figure(2)
for ii=1:(order+1)*2;
	subplot(2,order+1,ii)
	plot(theta(:,ii))
	title('Chain realisations') 
end
mean(theta)
theta_0'

theta_0'
