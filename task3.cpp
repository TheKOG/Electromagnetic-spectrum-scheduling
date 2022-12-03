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
int dfn;
bool Dfs(int id,double lim,int goal){
    //cerr<<id<<endl;
    if(goal==0)return true;
    if(goal>id)return false;
    if(id==n){
        dfn=0;
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
    dfn++;
    if(dfn>N)return false;
    for(int i=0;i<steps;i++,wave[id].phi+=step){
        //cerr<<"fkpps\n";
        double farthest=inf;
        for(int j=id+1;j<=n;j++){
            if(!wave[j].place)continue;
            double tmp=first_conflict(wave[id],wave[j]);
            farthest=min(farthest,tmp);
        }
        if(farthest<lim)continue;
        if(Dfs(id-1,lim,goal-1))return true;
    }
    wave[id].place=false;
    return Dfs(id-1,lim,goal);
}
Wave wave_[N];
bool change[N];
double prop;
double phi[N];
bool fuckpps=true;
bool Check(Wave ugt,double lim,double propmax,int boundary=100){
    for(int i=1;i<=n;i++)wave_[i]=wave[i],change[i]=false;
    int csum=0;
    for(int turn=1;turn<=boundary;turn++){
        int id=-1,maxc=0;
        for(int i=1;i<=n;i++){
            int conflicts=0;
            double fc=first_conflict(wave_[i],ugt);
            if(fc<lim)conflicts=n;
            for(int j=1;j<=n;j++){
                if(i==j)continue;
                fc=first_conflict(wave_[i],wave_[j]);
                if(fc>=lim)continue;
                conflicts++;
            }
            if(conflicts>maxc){
                maxc=conflicts;
                id=i;
            }
        }
        if(id==-1){
            prop=(double)csum/n;
            return prop<=propmax+eps;
        }
        if(!change[id]){
            csum++;
            change[id]=true;
        }
        if((double)csum/n>propmax)return false;
        //if(lim>368&&lim<369)cerr<<id;
        Wave pwave=wave_[id];
        for(int i=0;i<steps;i++){
            pwave.phi=pwave.T*i/steps;
            int conflicts=0;
            double fc=first_conflict(pwave,ugt);
            if(fc<lim)conflicts=n;
            for(int j=1;j<=n;j++){
                if(id==j)continue;
                fc=first_conflict(pwave,wave_[j]);
                if(fc>=lim)continue;
                conflicts++;
            }
            if(conflicts<maxc){
                maxc=conflicts;
                wave_[id]=pwave;
            }
        }
    }
    return false;
}
void Urgent_Task(double maxlim,double propmax=0.3){
    Wave ugt;
    scanf("%lf%lf",&ugt.D,&ugt.T);
    ugt.D-=ml;
    double l=0,r=maxlim,ansfc=0,ansphi=0,ansprop=0;
    while(r-l>minL*eps){
        double mid=(l+r)/2;
        bool flag=false;
        cerr<<"checking "<<mid<<"\n";
        for(int i=0;i<steps;i++){
            ugt.phi=ugt.T*i/steps;
            if(Check(ugt,mid,propmax)){
                if(prop<ansprop||flag==false){
                    ansphi=ugt.phi;
                    ansprop=prop;
                    for(int i=1;i<=n;i++){
                        while(wave_[i].phi>wave_[i].T)wave_[i].phi-=wave_[i].T;
                        phi[wave_[i].id]=wave_[i].phi;
                    }
                }
                flag=true;
            }
        }
        if(flag){
            cerr<<mid<<" is achievabel"<<endl;
            ansfc=l=mid;
        }else{
            cerr<<mid<<" is unachievabel"<<endl;
            r=mid;
        }
    }
    printf("With changed proportion of %.2lf the maximum first interference is %.2lf\n",ansprop,ansfc);
    for(int i=1;i<=n;i++){
        printf("phi(%d) is %.2lf\n",i,phi[i]);
    }
    printf("urgent phi is %.2lf\n",ansphi);
}
signed main(){
    //freopen("data/task3/1_1.txt","r",stdin);
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
    while(r-l>minL*eps*10){
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
        for(int j=1;j<=n;j++){
            if(wave[j].id==i)wave[j].phi=phi[i];
        }
    }
    printf("Now considering insert a urgent task.\n");
    Urgent_Task(ans);
    return 0;
}
/*
*/
