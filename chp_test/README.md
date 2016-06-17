
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **chp_test** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : chp_test

Published in : 'P. Burdejova, W.K. Härdle, P.Kokoszka and Q.Xiong (2015): Change point and trend
analyses of annual expectile curves of tropical storms'

Description : 'Change point test with hypothesis that mean functions of hurricane expectile curves
do not change with year'

Keywords : change point, test, curve, expectile, time varying, hurricane

See also : P_beta_est.R, data_load_hurricanes.R

Author : Burdejova P.

Submitted : 20160610 by Burdejova P.

Datafile : Hurricanes.csv

```


### R Code:
```r
# Run Data_load_hurricanes.R to obtain the object "expect_list"
# source (data_load_hurricanes.R)

#reaarange to have matrix of expectiles for given tau
expect_list_pertau = vector(mode='list', length=9)
exp_pertau	       = matrix(nrow=1460, ncol=65)

for (tau in 1:9){
	expect_list_pertau[[tau]] = exp_pertau
}

for (year in 1:65){
	year_data = expect_list[[year]]
	for (tau in 1:9){
		expect_list_pertau[[tau]][,year] = year_data[,tau]	
	}
}
		
# compute PC and scores
d=rep(0,9)
scores_list_pertau = vector(mode='list', length=9)
scores_pertau	   = matrix(nrow=1460, ncol=65)

for (tau in 1:9){
	scores_list_pertau[[tau]]= scores_pertau
}

eigenv_list_pertau = vector(mode='list', length=9)

for (tau in 1:9){	
	basis      = create.bspline.basis( c(0,1), nbasis=180, norder=4)
	fun_pec    = Data2fd(argvals=seq(0,1,length=1460),
			 		y=expect_list_pertau[[tau]], basisobj=basis)
	pc 	    	   = pca.fd(fun_pec, nharm=20)

	expl_var   = cumsum(pc$varprop)
	ind_over85 = which(expl_var > 0.85)
	d[tau]	   = min(ind_over85)
	d_tau      = min(ind_over85)
	
	pc_scores  				  = pc$scores
	expect_list_pertau[[tau]] = pc_scores	
	eigenv_list_pertau[[tau]] = pc$values	
}
			
# compute test statistic
S_tau = rep(0,9)
N	  = 65

for (tau in 1:9){	
	d_tau 	  = d[tau]
	pc_scores = expect_list_pertau[[tau]]	
	sum		  = 0
	for (l in 1:d_tau){
		sum_l=0
		for (k in 1:N){
			sum1  = sum(pc_scores[1:k,l])
			sum2  = sum(pc_scores[(k):N,l])
			diff  = sum1-( (k/N) * sum2)
			sum_l = sum_l + (diff*diff)
		} # of K
		sum_l = sum_l / (eigenv_list_pertau[[tau]][l])
		sum   = sum+sum_l
	} # of L
	S_tau[tau] = sum/(N*N)
} # of tau

print(d)
print(S_tau)
		
					
```
