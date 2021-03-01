#include <iostream>
#include <math.h>
#include <random>
#include <ctime>

using namespace std;

double gen_uniform(){
    std::random_device rd;
    std::mt19937_64 generator(rd());
    std::uniform_real_distribution<double> distribution;
    return distribution(generator);
}

double ksi(double w,double alpha){
    return pow(-log(w),0.25)/alpha;
}

double eta(double w){
    return sqrt(-log(w));
}

double gen_q_B(double a){
    //double w=gen_uniform();
    double k=ksi(gen_uniform(),a);
    double t=eta(gen_uniform());
    if(k<t){
        return 1;
    }
    else{
        return 0;
    }
}

double gen_q_C(double alpha){
    double k=ksi(gen_uniform(),alpha);
    return exp(-(k*k));
}

double gen_q_D(double alpha){
    double t=eta(gen_uniform());
    return 1-exp(-pow(alpha*t,4.0));
}

double gen_q_E(double alpha){
    double t1=-log(gen_uniform());
    double t2=-log(gen_uniform());
    double t3=-log(gen_uniform());
    double b = sqrt(t1+t2+t3);
    return 2*(1-exp(-pow(alpha*b,4.0)))/pow(b,4.0);
}

void task(double n_0,double alpha,double (*gen_q)(double)){
    clock_t begin=clock();
    double z2=pow(2.575,2.0);
    double e2=pow(0.01,2.0);
    double q=gen_q(alpha);
    double sum=q;
    double sum_square=q*q;
    double n=1;
    while (n<n_0){
        n+=1;
        q=gen_q(alpha);
        sum+=q;
        sum_square+=q*q;
    }
    if(sum==0.0){
        cout<<"Error: sum q is "<<sum<<".Need more bigger n_0";
        return;
    }
    double Q=sum/n;
    double sigma2=(sum_square-n*Q*Q)/(n-1);
    double criteria=(z2*sigma2)/(e2*Q*Q);
    while (n<criteria){
        n+=1;
        q=gen_q(alpha);
        sum+=q;
        sum_square+=q*q;
        Q=sum/n;
        sigma2=(sum_square-n*Q*Q)/(n-1);
        criteria=(z2*sigma2)/(e2*Q*Q);
    }
    clock_t end=clock();
    double estimated_seconds = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"alpha="<<alpha<<endl;
    cout<<"n* = "<<n<<endl;
    cout<<"Q = "<<Q<<endl;
    cout<<"Estimate time: "<<estimated_seconds<<" seconds"<<endl<<endl;
    
}

int main(int argc, char* argv[]){
    if(argc>1){
        double n_0=double(atoi(argv[2]));
        double alpha=double(atof(argv[3]));
        switch (atoi(argv[1])){
        case 1:
            task(n_0,alpha,gen_q_B);//n_0=10000;if alpha==0.1: n_0=1_000_000
            break;
        case 2:
            task(n_0,alpha,gen_q_C);//n_0=10000
            break;
        case 3:
            task(n_0,alpha,gen_q_D);//n_0=1000
            break;
        case 4:
            task(n_0,alpha,gen_q_E);//n_0=10000  
            break;
        }
    }
    else{
        cout<<"wrong args"<<endl;
    }
    
    
    return 0;
}