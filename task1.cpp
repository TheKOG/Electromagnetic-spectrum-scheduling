#include<bits/stdc++.h>
using namespace std;
const int N=114514;
const double eps=1e-3,inf=1145141919810;
int n,tot;
int steps=50;
struct Wave{
    double D,T,phi;
    int id;
    bool place;
    bool operator < (const Wave &a)const{return T<a.T;}
}wave[N];
double first_conflict(Wave a,Wave b,int boundary=100){
    if(b<a)swap(a,b);
    if(a.D<0||b.D<0)return inf;
    double delt=a.phi-b.phi,width=(a.D+b.D)/a.T,k=b.T/a.T;
    double U0=-delt/a.T+width/2,L0=-delt/a.T-width/2;
    double k1=k-floor(k),U0_=ceil(U0),L0_=floor(L0);
    if(k1<eps&&U0_-L0_<=1)return inf;
    if(k1<eps&&U0_-L0_>1+eps){
        double re=b.phi-floor(b.phi/b.T)*b.T;
        return re;
    }
    for(int i=-boundary;i<=boundary;i++){
        double F=(U0_-U0+i)/k1,G=(U0_-L0+i)/k1;
        int b_min=(int)ceil(F),b_max=(int)floor(G);
        if(b_max*b.T+b.phi<0)continue;
        for(int b_=b_min;b_<=b_max;b_++){
            double Ub=U0+k*b_;
            int a_=floor(Ub);
            if(b_*b.T+b.phi<0||a_*a.T+a.phi<0)continue;
            return b_*b.T+b.phi;
        }
    }
}
double maxL=-inf,minL=inf,ml;
bool Dfs(int id,double lim,int goal){
    //cerr<<id<<endl;
    if(goal==0)return true;
    if(goal>id)return false;
    if(id==n){
        wave[id].phi=0;
        wave[id].place=true;
        if(Dfs(id-1,lim,goal-1))return true;
        else{
            wave[id].place=false;
            return Dfs(id-1,lim,goal);
        }
    }
    wave[id].phi=wave[id+1].phi;
    wave[id].place=true;
    if(wave[id].D<0)
        return Dfs(id-1,lim,goal-1);
    double step=wave[id].T/steps;
    for(int i=0;i<steps;i++,wave[id].phi+=step){
        //cerr<<"fkpps\n";
        double farthest=inf;
        for(int j=id+1;j<=n;j++){
            if(!wave[j].place)continue;
            double tmp=first_conflict(wave[id],wave[j]);
            farthest=min(farthest,tmp);
            //printf("%d conflict %d is %lf   lim is %lf\n",id,j,tmp,lim);
        }
        if(farthest<lim)continue;
        if(Dfs(id-1,lim,goal-1))return true;
    }
    wave[id].place=false;
    return Dfs(id-1,lim,goal);
}
double phi[N];
signed main(){
    freopen("data/task1/2.txt","r",stdin);
    //freopen("result/task1/1.out","w",stdout);
    scanf("%d%lf",&n,&ml);
    double dsum=0;
    for(int i=1;i<=n;i++){
        scanf("%lf%lf",&wave[i].D,&wave[i].T);
        maxL=max(maxL,wave[i].T);
        wave[i].D-=ml;
        if(wave[i].D>0)minL=min(minL,wave[i].T);
        wave[i].id=i;
        dsum+=wave[i].D;
        //printf("%lf %lf\n",wave[i].D,wave[i].T);
    }
    sort(wave+1,wave+n+1);
    double l=0,r=maxL*128,ans=0;
    while(r-l>minL*eps){
        double mid=(r+l)/2;
        cerr<<"checking "<<mid<<"\n";
        if(Dfs(n,mid,n)){
            l=ans=mid;
            cerr<<mid<<" is achievabel"<<endl;
            for(int i=1;i<=n;i++){
                while(wave[i].phi>wave[i].T)wave[i].phi-=wave[i].T;
                phi[wave[i].id]=wave[i].phi;
                //printf("T(%d)=%lf\n",i,wave[i].T);
            }
        }else{
            r=mid;
            cerr<<mid<<" is unachievabel"<<endl;
        }
    }
    if(ans<eps){
        printf("There is no way to prevent the first interference from occuring in the first cycle.");
        return 0;
    }
    printf("The maximum time of the first interference is %.2lf\n",ans);
    for(int i=1;i<=n;i++){
        printf("phi(%d) is %.2lf\n",i,phi[i]);
    }
    return 0;
}
/*
0.04384768436222811 1.9675609671005667
0.020277296636061357 1.0164766671218217
*/
