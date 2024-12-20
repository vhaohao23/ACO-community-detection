#include<bits/stdc++.h>
#include <windows.h>
using namespace std;
const int N=500;//số thành phố tối đa
double p=0.8;
int A[N+1][N+1]={};
int numberCity;
const int maxInteration=100;
const int numberAnt=30  ;
int NE;
vector<vector<double>> pheromone ( N, vector<double>(N,1) );
double ALPHA=1,BETA=2;
double h[N+1][N+1]={};
double u[N]={};
double o[N]={};
double Tmax,Tmin;
int k[N]={};
int Community[N];
vector<int> e[N];
struct edge{
        int u,v;
};
double bestT,worstT;
vector<edge> trace[numberAnt+1];

int dpt=0;
int better=0;
vector<int> karate[100];
void solution(int k){
    double chance[N+1][N+1]={};
   for (int i=1;i<=numberCity;i++){
        double totalChance=0;
         //random để chọn thành phố tiếp theo
        for (int j=1;j<=numberCity;j++)
            if (i!=j&&A[i][j]){
                chance[i][j]=pow(pheromone[i][j],ALPHA)*pow(h[i][j],BETA);
                totalChance+=chance[i][j];
            }
       
        double randCity=(double(rand())/double(RAND_MAX))*totalChance;
        double check=0;
        int nextCity=0;
        for (int j=1;j<=numberCity;j++)
            if (i!=j&&A[i][j]){
                check+=chance[i][j];
                if (check>=randCity) {nextCity=j;break;}
            }
        
        e[i].push_back(nextCity);
        e[nextCity].push_back(i);
        trace[k].push_back({i,nextCity});
   }
}

void updatePheromone(double Modularity, int k,int iter){
    //pheromone bị bay hơi đi một lượng với tốc độ bay hơi là p

    for (int i=1;i<=numberCity;i++)
        for (int j=1;j<=numberCity;j++){    
            if (iter<=10)
            pheromone[i][j]*=p;
            
        }

    //cập nhật lại nồng độ pheromone trên đường khi con kiến đi qua
    
    bool dd[numberCity+1][numberCity+1]={};
    for (auto ve : trace[k]){
        int from=ve.u,to=ve.v;
            
        if (dd[from][to]) continue;
        pheromone[from][to]+=Modularity;
        pheromone[to][from]+=Modularity;
            
        dd[from][to]=dd[to][from]=1;                                                                                                 
    }

    for (int i=1;i<=numberCity;i++)
        for (int j=1;j<=numberCity;j++){
            bestT=max(bestT,pheromone[i][j]),worstT=min(worstT,pheromone[i][j]);
            if (pheromone[i][j]>Tmax) {pheromone[i][j]=pheromone[j][i]=Tmax;better++;}
            if (pheromone[i][j]<Tmin) pheromone[i][j]=pheromone[j][i]=Tmin;

        }
}

void bfs(){
    int cnt=0;
    int dd[numberCity+1]={};
    for (int i=1;i<=numberCity;i++)
        if (!dd[i]){
            ++cnt;
            queue<int> q;
            q.push(i);
            dd[i]=1;
            while (!q.empty())
            {
                int u=q.front();
                q.pop();
                Community[u]=cnt;
                for (int v:e[u])
                if (!dd[v]){
                    dd[v]=1;
                    q.push(v);
                }
            }

        }
}

double measureCrit(){
    double q=0;
    for (int i=1;i<=numberCity;i++)
        for (int j=1;j<=numberCity;j++)
        // if (i!=j)
            q+=((double(A[i][j]) - double(k[i]*k[j])/double(2.0*NE) )*(Community[i]==Community[j]));

    return q/double(2.0*NE);
}

void setting(){
    double ans=-1;
    int ansCommunity[numberCity+1]={};
    for (int i=1;i<=maxInteration;i++){
        bestT=0,worstT=100;

        double modularity[numberAnt]={};
        double ib=-1;
        int whichAnt;

        for (int k=0;k<numberAnt;k++){
            //rs
            for (int j=1;j<=numberCity;j++) e[j].clear();

            solution(k);

            memset(Community,sizeof(Community),0);
            bfs();

            modularity[k]=measureCrit();

            if (modularity[k]>ans){
                for (int node=1;node<=numberCity;node++)
                    ansCommunity[node]=Community[node];

                ans=modularity[k];

            }

            if (modularity[k]>ib) ib=modularity[k],whichAnt=k;

        }

        Tmax=ans/(1-p),Tmin=Tmax*0.005;
        // cout<<trace[whichAnt].size()<<" "<<ib<<"\n";
        updatePheromone(ib,whichAnt,i);

        for (int k=0;k<numberAnt;k++) trace[k].clear();

        cout<<Tmin<<" "<<ib<<" "<<Tmax<<" "<<bestT<<","<<worstT<<"\n";
        // cout<<better<<"\n";

    }
    cout<<"Các nhóm:\n";
    for (int i=1;i<=numberCity;i++)
        karate[ansCommunity[i]].push_back(i);

    for (auto x:karate){
        for (auto y:x) cout<<y<<" ";
        if (x.size()>0)
        cout<<"\n";
    }
    cout<<"Modularity global best:";
    cout<<ans<<"\n";

    // cout<<"Membership vector:\n";
    // for (int i=1;i<=numberCity;i++)
    //     cout<<ansCommunity[i]<<" ";


    //  for (int i=1;i<=numberCity;i++){
    //     for (int j=1;j<=numberCity;j++)
    //         cout<<pheromone[i][j]<<" ";
    //     cout<<"\n";
    // }
        
    // cout<<dpt;
    // cout<<better;
}

//calculate heruristic in4
void calU(int i){
    for (int l=1;l<=numberCity;l++)
        u[i]+=(double(A[i][l])/double(numberCity));
}
void calO(int i){
    for (int l=1;l<=numberCity;l++)
        o[i]+=(pow(A[i][l]-u[i],2)/double(numberCity));


    o[i]=sqrt(o[i]);
}
double C(int i,int j){
    double sum=0;
    for (int l=1;l<=numberCity;l++)
    {
        sum+=( (double(A[i][l])-u[i])*(double(A[j][l])-u[j]) );
    }
    sum/=(double(numberCity)*o[i]*o[j]);
    return sum;
}
void calH(int i,int j){
    h[i][j]=h[j][i]=1/(1+exp(-C(i,j)));
}
void heuristic(){
    for (int i=1;i<=numberCity;i++)
        calU(i);
    for (int i=1;i<=numberCity;i++)
        calO(i);

    //cal heuristic in4
    for (int i=1;i<=numberCity;i++)
        for (int j=1;j<=numberCity;j++)
            if (i!=j)
                calH(i,j);

}

int main(){
    freopen("population.txt","r",stdin);
     srand((unsigned int)time(NULL));

    cin>>numberCity;
    cin>>NE;
    for (int i=1;i<=NE;i++){
        int u,v,w;
        cin>>u>>v;
        
        w=1;
        if (!A[u][v])
            A[u][v]=A[v][u]=w;

    }
    //calculate degree
    for (int i=1;i<=numberCity;i++)
    {

        for (int j=1;j<=numberCity;j++)
            k[i]+=A[i][j];
        // A[i][i]=1;
    }
    heuristic();
    
    setting();

}