// Medcouple from phylter / mrfDepth (July 2023).
// P. Segaert, based on M. Maechler and G. Brys. GPL-2-or-later.
#include "phylter2/statistics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
namespace phylter2 {
namespace {
using namespace std;
double sign(double x){
	return((x>0.0)?1.0:((x==0.0)?0.0:-1.0));
}
double h_kern(double a,double b,int ai,int bi,int ab,double eps){
    	if(fabs(a-b)<2.0*eps || b>0)	return sign((double)(ab-(ai+bi)));
    	return(a+b)/(a-b);
}
double whimed_i(double a[],int w[],int n,double a_cand[],double a_srt[],int w_cand[]){
	int n2,i,kcand; /* sum of weights: `int' do overflow when  n ~>= 1e5 */
	int64_t wleft,wmid,wright,w_tot=0,wrest=0;
	double trial;

	for(i=0;i<n;++i)	w_tot+=w[i];
  	do{/* REPEAT : */
		for(i=0;i<n;++i)	a_srt[i]=a[i];
		n2=n/2;/*=^= n/2 +1 with 0-indexing */
		std::nth_element(a_srt,a_srt+n2,a_srt+n);//rPsort(a_srt,n,n2);
		trial=a_srt[n2];
		wleft=0;
		wmid=0;
		wright=0;
		for(i=0;i<n;++i){
		    	if(a[i]<trial)		wleft+=w[i];
		    	else if(a[i]>trial)	wright+=w[i];
		    	else			wmid+=w[i];
		}
		kcand=0;
		if(2*(wrest+wleft)>w_tot){
			for(i=0;i<n;++i){
				if(a[i]<trial){
					a_cand[kcand]=a[i];
				    	w_cand[kcand]=w[i];	
					++kcand;
				}
		    	}
		}
		else if(2*(wrest+wleft+wmid)<=w_tot){
			for(i=0;i<n;++i){
				if(a[i]>trial){
				    	a_cand[kcand]=a[i];
				    	w_cand[kcand]=w[i];
				    	++kcand;
				}
		    	}
		    	wrest+=wleft+wmid;
		} else {
		    	return trial;
		}
		n=kcand;
		for(i=0;i<n;++i){
			a[i]=a_cand[i];
		    	w[i]=w_cand[i];
		}
    	} while(1);
}
double mlmccN2(double x[],const int n){
	double medc,xmed,xden,x_eps,trial=0.0,eps_trial;
	//double eps[2],work[n];
	double eps[2];
	std::vector<double> work_storage(n);
	double * work = work_storage.data();
    int i,j,k,it=0,ind=n/2,h1,h2;
	//int iwt[n];
	std::vector<int> iwt_storage(n);
	int *iwt = iwt_storage.data();
	int64_t knew,nr,sum_p,sum_q,nl=0,neq=0;
    	bool converged=true,IsFound=false;

    	
	eps[0]=2.220446e-16;
	eps[1]=2.225074e-308;
    	if(n<3)	return(0.0);
    	if(n%2){ /* n even */
		xmed=x[ind+1];
    	} else { /* n  odd */
		xmed=(x[ind]+x[ind+1])/2;
    	}
    	if(fabs(x[1]-xmed)<eps[0]*(eps[0]+fabs(xmed))){
		return(-1.0);
    	} else if(fabs(x[n]-xmed)<eps[0]*(eps[0]+fabs(xmed))){
		return(1.0);
    	}
    	for(i=1;i<=n;i++)	x[i]-=xmed;
    	xden=-1.0/(2.0*max(-x[1],x[n]));
    	for(i=1;i<=n;i++)	x[i]*=xden;
    	xmed*=xden;
    	j=1;
    	x_eps=eps[0]*(eps[0]+fabs(xmed));
	while(j<=n && x[j]>x_eps)	j++;
    	i=1;
    	double *x2=x+j-1;/* pointer -- corresponding to  x2[i]=x[j];*/
	while(j<=n && x[j]>-x_eps){ /* test relative to xmed */
    	    	j++;
    	    	i++;
    	}
    	h1=j-1;/*== size of x1[]== the sum of those two sizes above */
    	h2=i+(n-j);//== size of x2[]== maximal size of whimed()arrays

	double * acand;
    double * a_srt;
    std::vector<double> acand_storage(h2);
    acand = acand_storage.data();
    std::vector<double> a_srt_storage(h2);
    a_srt = a_srt_storage.data();
    int * iw_cand;
    int * left;
    int * right;
    int * p; 
    int * q;
    std::vector<int> iw_cand_storage(h2);
    iw_cand = iw_cand_storage.data();
    std::vector<int> left_storage(h2+1);
    left = left_storage.data();
    std::vector<int> right_storage(h2+1);
    right = right_storage.data();
    std::vector<int> p_storage(h2+1);
    p = p_storage.data();
    std::vector<int> q_storage(h2+1);
    q = q_storage.data();

	std::fill_n(left,h2+1,1);
	std::fill_n(right,h2+1,h1);
    	nr=((int64_t)h1)*((int64_t)h2);// <-- careful to *NOT* overflow
    	knew=nr/2+1;
    	while(!IsFound && (nr-nl+neq>n) && it<1000){
		it++;
		j=0;
		for(i=1;i<=h2;i++){
		    	if(left[i]<=right[i]){
				iwt[j]=right[i]-left[i]+1;
				k=left[i]+(iwt[j]/2);
				work[j]=h_kern(x[k],x2[i],k,i,h1+1,eps[1]);
				j++;
		    	}
		}
		trial=whimed_i(work,iwt,j,acand,a_srt,iw_cand);
		eps_trial=eps[0]*(eps[0]+fabs(trial));
		j=1;
		for(i=h2;i>=1;i--){
		    	while(j<=h1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1])-trial>eps_trial)	j++;
		    	p[i]=j-1;
		}
		j=h1;
		for(i=1,sum_p=0,sum_q=0;i<=h2;i++){
		    	while(j>=1 && trial-h_kern(x[j],x2[i],j,i,h1+1,eps[1])>eps_trial)	j--;
		    	q[i]=j+1;
		    	sum_p+=p[i];
		    	sum_q+=j;/*=q[i]-1 */
		}
		if(knew<=sum_p){
		    	for(i=1,neq=0;i<=h2;i++){
				right[i]=p[i];
				if(left[i]>right[i]+1)neq+=left[i]-right[i]-1;
		    	}
		    	nr=sum_p;
		} else { /* knew>sum_p */
		    	IsFound=(knew<=sum_q);/* i.e. sum_p < knew<=sum_q */;
		    	if(IsFound){
			 	medc=trial;
		    	} else { /*	 knew>sum_q */
		    	    	for(i=1;i<=h2;i++){
				    	left[i]=q[i];
				    	if(left[i]>right[i]+1)neq+=left[i]-right[i]-1;
				}
				nl=sum_q;
		    	}
		}
    	} /* end while loop */
    	converged=IsFound || (nr-nl+neq<=n);
    	if(!converged)	medc=trial;
    	if(converged && !IsFound){ /* e.g.,for  mc(1:4): */
		j=0;
		for(i=1;i<=h2;i++){
		    	if(left[i]<=right[i]){
				for(k=left[i];k<=right[i];k++){
				    	work[j]=-h_kern(x[k],x2[i],k,i,h1+1,eps[1]);
				    	j++;
				}
		    	}
		}
		std::nth_element(work,work+knew-nl-1,work+j);
		medc=-work[knew-nl-1];
    	}

        
        
        
        
        
        
        
                

    	return medc;
} 
double mlmccN(const double z[],const int n,const int dr){
	//double x[n+1];
    double * x;
    std::vector<double> x_storage(n+1);
    x = x_storage.data();
	double res,res_dr;
	int i;

    x[0]=0.0;
    for(i=0;i<n;i++)	x[i+1]=-z[i];
	std::sort(x+1,x+n+1);
	if(dr){
		//double w[n+1];
        double * w;
        std::vector<double> w_storage(n+1);
    w = w_storage.data();
	    	w[0]=0.0;
	    	for(i=0;i<n;i++)	w[i+1]=-x[n-i]-x[1];
		res_dr=mlmccN2(w,n);
	 	
	}
	res=mlmccN2(x,n);
	if(dr)	res=(res-res_dr)/2.0;
    
	return(res);
}

} // namespace
double medcouple(const std::vector<double>& values) {
  if (values.empty() || values.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("unsupported medcouple sample size");
  for (double value : values) if (!std::isfinite(value))
    throw std::invalid_argument("medcouple requires finite values");
  return mlmccN(values.data(), static_cast<int>(values.size()), values.size() > 100);
}
} // namespace phylter2

