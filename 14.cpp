#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

double expSum()
{
	int i=0;
	double sum=0;
	while (i<12)
	{
		sum+=rand()*1.0/RAND_MAX;
		i++;
	}
	return sum-6;
}

double precSum()
{
	double s=expSum();
	return s+(s*s*s+3*s)/240;
}

double BoxMuller()
{
	double s=sqrt(-2*log(rand()*1.0/RAND_MAX))*cos(2*3.14*rand()*1.0/RAND_MAX);
	return s;
}

double av(double *probs,int n)
{
	double avrg=0;
	int i=0;
	while(i<n)
	{
		avrg+=i*probs[i];
		i++;
	}
	return avrg;
}


double var(double *probs,int n)
{
	double variance=0;
	int i=0;
	while(i<n)
	{
		variance+=i*i*probs[i];
		i++;
	}
	variance=variance-av(probs,n)*av(probs,n);
	return variance;
}

int main()
{
	srand(unsigned(time(0)));
	double avr,variance,num;
	cout<<"Enter average,variance and number of experiments:"<<endl;
	cin>>avr;
	cin>>variance;
	cin>>num;
	int n=int(ceil(sqrt(num)));
	double *theoretic=new double[n];
	double intStart=avr-3*variance;
	double intLen=6*variance/n;
	int i=0;
	double pr=1;
	while (i<n)
	{
		double center=intStart+intLen*(i+0.5);
		theoretic[i]=intLen*exp(-(center-avr)*(center-avr)/(2*variance*variance))/(variance*sqrt(2*3.14159));
		pr-=theoretic[i];
		i++;
	}
	theoretic[i-1]+=pr/2;
	theoretic[0]+=pr/2;
	double *expProb=new double[n];
	cout<<"Sum method:"<<endl;
		i=0;
		while (i<n)
		{
			expProb[i]=0;
			i++;
		}
		i=0;
		while (i<num)
		{
			double result=variance*expSum()+avr;
			result-=intStart;
			result/=intLen;
			if (result<0)
			   result=0;
	        if (result>=n)
	           result=n-0.5;
			expProb[int(floor(result))]+=1.0/num;
			i++;
		}
		cout<<endl<<"Experemental probabilities:"<<endl;	
		i=0;
		while (i<n)
		{
			cout<<intStart+i*intLen<<" to "<<intStart+(i+1)*intLen<<": "<<expProb[i]<<endl;
			i++;
		}
		double error=abs(av(expProb,n)-av(theoretic,n))/(av(theoretic,n)+1);
		cout<<"Average: "<<av(expProb,n)+1<<" ( "<<error*100<<"% error)"<<endl;
		
		error=abs(var(expProb,n)-var(theoretic,n))/var(theoretic,n);
		cout<<"Variance: "<<var(expProb,n)+1<<" ( "<<error*100<<"% error)"<<endl;
		
		double chisqr=0;
		i=0;
		while (i<n)
		{
			chisqr+=(expProb[i]*expProb[i])/theoretic[i];
			i++;
		}
		chisqr=(chisqr-1)*num;
		cout<<"Chi squared: "<<chisqr<<endl<<"Compare this to table value to test destribution"<<endl;
		cout<<"(Smaller means that experemental matches theoretical)";
	cout<<endl<<endl<<"Precise sum method:"<<endl;
		i=0;
		while (i<n)
		{
			expProb[i]=0;
			i++;
		}
		i=0;
		while (i<num)
		{
			double result=variance*precSum()+avr;
			result-=intStart;
			result/=intLen;
			if (result<0)
			   result=0;
	        if (result>=n)
	           result=n-0.5;
			expProb[int(floor(result))]+=1.0/num;
			i++;
		}
		cout<<endl<<"Experemental probabilities:"<<endl;	
		i=0;
		while (i<n)
		{
			cout<<intStart+i*intLen<<" to "<<intStart+(i+1)*intLen<<": "<<expProb[i]<<endl;
			i++;
		}
		error=abs(av(expProb,n)-av(theoretic,n))/(av(theoretic,n)+1);
		cout<<"Average: "<<av(expProb,n)+1<<" ( "<<error*100<<"% error)"<<endl;
		
		error=abs(var(expProb,n)-var(theoretic,n))/var(theoretic,n);
		cout<<"Variance: "<<var(expProb,n)+1<<" ( "<<error*100<<"% error)"<<endl;
		
		chisqr=0;
		i=0;
		while (i<n)
		{
			chisqr+=(expProb[i]*expProb[i])/theoretic[i];
			i++;
		}
		chisqr=(chisqr-1)*num;
		cout<<"Chi squared: "<<chisqr<<endl<<"Compare this to table value to test destribution"<<endl;
		cout<<"(Smaller means that experemental matches theoretical)";
	cout<<endl<<endl<<"Box-Muller method:"<<endl;
		i=0;
		while (i<n)
		{
			expProb[i]=0;
			i++;
		}
		i=0;
		while (i<num)
		{
			double result=variance*BoxMuller()+avr;
			result-=intStart;
			result/=intLen;
			if (result<0)
			   result=0;
	        if (result>=n)
	           result=n-0.5;
			expProb[int(floor(result))]+=1.0/num;
			i++;
		}
		cout<<endl<<"Experemental probabilities:"<<endl;	
		i=0;
		while (i<n)
		{
			cout<<intStart+i*intLen<<" to "<<intStart+(i+1)*intLen<<": "<<expProb[i]<<endl;
			i++;
		}
		error=abs(av(expProb,n)-av(theoretic,n))/(av(theoretic,n)+1);
		cout<<"Average: "<<av(expProb,n)+1<<" ( "<<error*100<<"% error)"<<endl;
		
		error=abs(var(expProb,n)-var(theoretic,n))/var(theoretic,n);
		cout<<"Variance: "<<var(expProb,n)+1<<" ( "<<error*100<<"% error)"<<endl;
		
		chisqr=0;
		i=0;
		while (i<n)
		{
			chisqr+=(expProb[i]*expProb[i])/theoretic[i];
			i++;
		}
		chisqr=(chisqr-1)*num;
		cout<<"Chi squared: "<<chisqr<<endl<<"Compare this to table value to test destribution"<<endl;
		cout<<"(Smaller means that experemental matches theoretical)";

}
