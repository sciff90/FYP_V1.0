function pvalue = ptarget(y,u,elim,theta,order)
	
	yp = filter(theta((order+2):2*(order+1)),theta(1:order+1),u);
	penew = y-yp;
	
	p1 = max(abs(penew)>elim);
	
	  if p1>0  % At least one residual violated uniform bound
            pvalue = -1e100;  % Close enough to minus infinity I guess
      else
            pvalue = -length(penew)*log(elim);
      end;
	
end
