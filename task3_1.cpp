#include<bits/stdc++.h>
using namespace std;
const double eps=1e-5;
const int steps=100;
struct spec{
    double phi,T,d,hz;
}a[111];
int n,cnt,rid;
double len,ml,mt=0,md=0,Hz[16],rans=0,reward=0,rT,K=0.04,K2=0.001;
inline bool cmp(spec v1,spec v2){
    return v1.hz<v2.hz;
}
vector<int> p[16];
vector<spec> ans[101];
struct F{
    int lstp=0,lsts=0;
    double du=0,phi=0,val=0;
}f[steps+20][1<<15];
double first_conflict(spec A,spec b,int boundary=100){
    if(b.T<A.T) swap(A,b);
    if(A.d<0||b.d<0) return -1;
    double delt=A.phi-b.phi,width=(A.d+b.d)/A.T,k=b.T/A.T;
    double U0=-delt/A.T+width/2,L0=-delt/A.T-width/2;
    double k1=k-floor(k),U0_=ceil(U0),L0_=floor(L0);
    if(k1<eps&&U0_-L0_<=1)return -1;
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
            if(b_*b.T+b.phi<0||a_*A.T+A.phi<0)continue;
            return b_*b.T+b.phi;
        }
    }
}

double count_conflict(spec A,spec b,int boundary=100){
    if(b.T<A.T) swap(A,b);
    if(A.d<0||b.d<0) return 0;
    double delt=A.phi-b.phi,width=(A.d+b.d)/A.T,k=b.T/A.T;
    double U0=-delt/A.T+width/2,L0=-delt/A.T-width/2;
    double k1=k-floor(k),U0_=ceil(U0),L0_=floor(L0);
    if(k1<eps&&U0_-L0_<=1)return 0;
    double re=0;
    if(k1<eps&&U0_-L0_>1+eps){
        re=floor((len-b.phi)/b.T);
        return re;
    }
    double first=0;
    for(int i=-boundary;i<=boundary;i++){
        double F=(U0_-U0+i)/k1,G=(U0_-L0+i)/k1;
        int b_min=(int)ceil(F),b_max=(int)floor(G);
        if(b_max*b.T+b.phi<0)continue;
        if(b_min*b.T+b.phi>len)break;
        for(int b_=b_min;b_<=b_max;b_++){
            double Ub=U0+k*b_;
            int a_=floor(Ub);
            if(b_*b.T+b.phi<0||a_*A.T+A.phi<0)continue;
            if(b_*b.T+b.phi>len||a_*A.T+A.phi>len)break;
            re+=1;
            if(first<eps)first=b_*b.T+b.phi;
        }
    }
    if(re>eps)re+=reward;
    return re*(2*len-first)/(2*len);
}
double conflict(double pi,double t,int id,double d){
    spec x;
    x.phi=pi;x.T=t;x.d=d;
    double re=count_conflict(x,a[id]);
    return re;
}
int count(int val){
    int re=0;
    while(val){
        ++re;
        val=val>>1;
    }
    return re;
}
signed main(){
    scanf("%d%lf%lf",&n,&len,&ml);
    for(int i=1;i<=n;i++){
        scanf("%lf%lf%lf%lf",&a[i].phi,&a[i].d,&a[i].T,&a[i].hz);
        mt=max(mt,a[i].T);
        md=max(md,a[i].d);
    }
    sort(a+1,a+n+1,cmp);
    cnt=0;
    bool flag=false;
    for(int i=1;i<=n;i++){
        a[i].d-=ml;
        if(a[i].d<0)continue;
        if(!flag){
            flag=true;
            p[++cnt].push_back(i);
        }else if(abs(a[i].hz-Hz[cnt])<eps){
            p[cnt].push_back(i);
        }else{
            p[++cnt].push_back(i);
        }
        Hz[cnt]=a[i].hz;
    }
    for(double t=0.1;t<=3;t+=0.1){
        double T=t*mt;
        for(int pv=1;pv<=steps;pv++){
            for(int j=0;j<(1<<cnt);j++){
                for(int i=1;i<pv;i++){
                    double d=(double)(pv-i)/steps*T;
                    double phi_now=(double)(pv+i)/(2*steps)*T;
                    if(f[pv][j].val<=f[i][j].val){
                        f[pv][j].val=f[i][j].val;
                        f[pv][j].du=d;
                        f[pv][j].phi=phi_now;
                        f[pv][j].lstp=i;
                        f[pv][j].lsts=j;
                    }
                    if(d<ml||d>md)continue;
                    for(int k=1;k<=cnt;k++){
                        double re=0;
                        if(!(j&(1<<(k-1))))continue;
                        for(auto z:p[k])
                            re+=conflict(phi_now,t*mt,z,d-ml);
                        int state=j^(1<<(k-1));
                        if(f[pv][j].val<=f[i][state].val+re+reward){
                            f[pv][j].val=f[i][state].val+re+reward;
                            f[pv][j].du=d;
                            f[pv][j].phi=phi_now;
                            f[pv][j].lstp=i;
                            f[pv][j].lsts=state;
                        }
                    }
                }
            }
        }
        int rp,rs;
        double re=0;
        for(int i=1;i<=steps;i++)
            for(int j=0;j<(1<<cnt);j++)
                if(re<f[i][j].val){
                    rp=i,rs=j;
                    re=f[i][j].val;
                }
        if(re<eps) continue;
        int tid=floor(0.1+t/0.1);
        while(1){
            int v1=f[rp][rs].lstp,v2=f[rp][rs].lsts;
            spec tmp;
            tmp.hz=Hz[count(rs^v2)];
            tmp.phi=f[rp][rs].phi;
            tmp.d=f[rp][rs].du;
            tmp.T=t*mt;
            if(rs!=v2){
                ans[tid].push_back(tmp);
            }
            rp=v1,rs=v2;
            if(!v2) break;
        }
        if(rans<re){
            rans=re;
            rid=tid;
            rT=t*mt;
        }
        memset(f,0,sizeof(0));
    }
    int m;
    scanf("%d",&m);
    for(int i=1;i<=m;i++){
        spec x;
        scanf("%lf%lf%lf%lf",&x.phi,&x.d,&x.T,&x.hz);
        flag=0;
        if(x.d<ml){
            printf("It's impossible to detect the signal\n");
            continue;
        }
        spec y;
        y.d=ml;y.T=rT;y.hz=x.hz;
        y.phi=x.phi-x.d/2+ml/2;
        double l=y.phi,r=x.phi+x.d/2-ml/2,mw=0,mphi=-1;
        for(auto v:ans[rid]){
            if(abs(v.hz-x.hz)>eps) continue;
            double dd;
            if(x.T<rT) dd=v.d;
            else dd=x.d;
            v.d-=ml/2;x.d-=ml/2;
            double fc=first_conflict(v,x);
            if(fc>eps) mw=K*(1-(fc-dd/2)/len)+K2;
        }
        for(double now_phi=l;now_phi<=r;now_phi+=0.01*(r-l)){
            double w=0,re=0;
            w+=(1-(now_phi-ml/2)/len)*K;
            y.phi=now_phi;
            if(y.phi>=rT){
                int kl=y.phi/rT;
                y.phi=y.phi-kl*rT;
            }
            for(auto v:ans[rid])
                if(count_conflict(v,y)) re=re+1;
            w+=K2*(1-re/ans[rid].size());
            if(mw<w) mw=w,mphi=now_phi;
        }
        if(mphi==-1) printf("Value:%lf\n No Change Required\n",mw);
        else{
            y.phi=mphi;
            for(int i=0;i<ans[rid].size();i++)
                if(count_conflict(ans[rid][i],y))
                    ans[rid][i].d=abs(y.phi-ml/2-ans[rid][i].phi);
            ans[rid].push_back(y);
            printf("We make changes like these:\n");
            printf("Value:%lf\n",mw);
            for(auto v:ans[rid])
                printf("phi=%lf duration=%lf frequency=%lf T=%lf\n",v.phi,v.d,v.hz,rT);
        }
    }
    return 0;
}
