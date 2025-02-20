#include<bits/stdc++.h>
using namespace std;
const int T=100;
int pop=100;
int N;
int NE;
double ALPHA=1,BETA=2;
double ans=0;
double pheMax,pheMin;
vector<vector<int>> A;
vector<vector<double>> h;
vector<vector<int>> e;
vector<double> u;
vector<double> o;
vector<int> k;
vector<int> bestCom;
vector<vector<double>> pheromone;
vector<vector<int>> trace;
double evaRate=0.8;
double epsilion;
random_device rd;   
mt19937 gen(rd());
vector<bool> dd;
vector<vector<int>> g(pop+1);
vector<vector<int>> dk(pop+1);
vector<vector<int>> lk(pop+1);
vector<vector<int>> x(pop+1);

void calU(int i){
    for (int l=1;l<=N;l++)
        u[i]+=(double(A[i][l])/double(N));
}
void calO(int i){
    for (int l=1;l<=N;l++)
        o[i]+=(pow(A[i][l]-u[i],2)/double(N));


    o[i]=sqrt(o[i]);
}
double C(int i,int j){
    double sum=0;
    for (int l=1;l<=N;l++)
    {
        sum+=( (double(A[i][l])-u[i])*(double(A[j][l])-u[j]) );
    }
    sum/=(double(N)*o[i]*o[j]);
    return sum;
}

void calH(int i,int j){
    h[i][j]=h[j][i]=1/(1+exp(-C(i,j)));
}
void heuristic(){
    for (int i=1;i<=N;i++)
        calU(i);
    for (int i=1;i<=N;i++)
        calO(i);

    for (int i=1;i<=N;i++)
        for (int j:e[i])
                calH(i,j);
}
void solution(vector<int> &select){
    double totalChance,randChance;
    for (int i=1;i<=N;i++)
        if (select[i]==i)
            {
                totalChance=0;
                for (int j:e[i])
                    totalChance+=pow(pheromone[i][j],ALPHA)*pow(h[i][j],BETA);
                
                uniform_real_distribution dis(0.0,totalChance);
                randChance=dis(gen);
                double check=0;

                for (int j:e[i]){
                    check+=pow(pheromone[i][j],ALPHA)*pow(h[i][j],BETA);
                    if (check>randChance) {select[i]=j;break;}
                }
            }
}
void decoding(vector<int> select,vector<int> &cs){
    trace.clear();
    trace.resize(N+1);
    
    for (int i=1;i<=N;i++)
        {
            trace[i].push_back(select[i]);
            trace[select[i]].push_back(i);
        }

    
    queue<int> q;
    int cnt=0;
    
    for (int i=1;i<=N;i++) dd[i]=0;

    for (int i=1;i<N;i++)
        if (!dd[i]){
            ++cnt;
            q.push(i);
            dd[i]=1;
            while (!q.empty()){
                int u=q.front();
                q.pop();
                cs[u]=cnt;
                for (int v:trace[u])
                    if (!dd[v]){
                        dd[v]=1;
                        q.push(v);
                    }
            }
        }
}
double modularity(vector<int> l){
    double Q=0;

    for (int i=1;i<=N;i++)
        for (int j=1;j<=N;j++)
                Q+=double(A[i][j]-k[i]*k[j]/double(2*NE))*double(l[i]==l[j]);
    return Q/double(2*NE);
}
void updatePheromone(vector<int> select,double qibest){
    for (int i=1;i<=N;i++)
        for (int j:e[i])
            pheromone[i][j]*=evaRate;
    int j;
    for (int i=1;i<=N;i++){
        j=select[i];
        pheromone[i][j]+=qibest;
        if (pheromone[i][j]>pheMax) pheromone[i][j]=pheMax;
        if (pheromone[i][j]<pheMin) pheromone[i][j]=pheMin;
    }
}
// void EPD(){
//     if (g.size()<10) return;

//     vector<pair<double, int>> modularityValues;
//     for (int i = 1; i <= pop; i++) {
//         double modValue = modularity();  
//         modularityValues.push_back({modValue, i});
//     }

//     sort(modularityValues.begin(), modularityValues.end());

//     vector<vector<int>> sortedG(pop + 1);
//     vector<vector<int>> sorteddk(pop + 1);
//     vector<vector<int>> sortedlk(pop + 1); 
//     for (int i = 0; i < pop; i++) {
//         sortedG[i + 1] = g[modularityValues[i].second];
//         sorteddk[i + 1] = dk[modularityValues[i].second];
//         sortedlk[i + 1] = lk[modularityValues[i].second];
//     }

//     g=sortedG;
//     dk=sorteddk;
//     lk=sortedlk;

//     double N_nor=pop-(pop/2+1)+1;
//     uniform_real_distribution<double> dis(0,1);
//     for (int i=pop/2+1;i<=pop;i++){
//         double C=1.0-exp(-double(i)/N_nor);
//         double rand=dis(gen);
//         if (rand<=C){
//             g.erase(g.begin() + i);
//             dk.erase(dk.begin() + i);
//             lk.erase(lk.begin() + i);   
//             --pop;
//         }
//     }
    
// }
void intersection(vector<int> a,vector<int> b,vector<int> &inter){
    for (int i=1;i<=N;i++)
        if (a[i]==b[i]) inter[i]=a[i];
}

void selfLearning(int p){
    vector<int> tmp;
    for (int i=1;i<=N;i++)
        for (int j:e[i]){
            tmp=g[p];
            g[p][i]=g[p][j];

            if (modularity(g[p])<modularity(tmp)){
                g[p]=tmp;
            }
        }
}

void ACO(){
    vector<pair<int,int>> save;
    vector<int> select;
    select.resize(N+1);
    vector<int> cs;
    cs.resize(N+1);
    vector<int> ibest;
    vector<int> selectib;
    vector<int> inter;
    inter.resize(N+1,0);
    int posOptimal=0,posSuboptimal=0;
    for (int t=1;t<=T;t++){
        double qibest=0;
       
        for (int p=1;p<=pop;p++){
            for (int i=1;i<=N;i++) select[i]=i,cs[i]=0;
            
            if (p>40)
            solution(select);
            else {
                solution(inter);
                select=inter;
            }
            x[p]=select;

            decoding(select,cs);
            g[p]=cs;
        }

        for (int p=1;p<=pop;p++){
            double q=modularity(g[p]);
            save.push_back({q,p});
        }
        
        sort(save.rbegin(), save.rend());
        posOptimal=save[0].second;
        posSuboptimal=save[1].second;
        selfLearning(posOptimal);
        selfLearning(posSuboptimal);

        for (int p=1;p<=pop;p++){
            double q=modularity(g[p]);

            if (q>qibest){
                qibest=q;
                ibest=g[p];
                selectib=x[p];
            }
        }

        if (qibest>ans){
            ans=qibest;
            bestCom=ibest;
        }

        for (int i=1;i<=N;i++) inter[i]=i;
        intersection(x[posOptimal],x[posSuboptimal],inter);

        pheMax=ans/(1-evaRate),pheMin=pheMax*epsilion;
        updatePheromone(selectib,qibest);
        cout<<ans<<"\n";
        // EPD();
    }

    cout<<ans<<"\n";
    for (int l:bestCom) cout<<l<<" ";
}
int main(){
    clock_t tStart = clock();

    freopen("input.txt","r",stdin);
    cin>>N>>NE;
    A.resize(N+1);
    for (int i=1;i<=N;i++) A[i].resize(N+1,0);
    u.resize(N+1,0);
    o.resize(N+1,0);
    k.resize(N+1,0);
    e.resize(N+1);
    h.resize(N+1);
    for (int i=1;i<=N;i++) h[i].resize(N+1,0);
    pheromone.resize(N+1);
    for (int i=1;i<=N;i++)
        pheromone[i].resize(N+1,1);
    trace.resize(N+1);
    dd.resize(N+1);
    for (int i=1;i<=pop;i++) {
        g[i].resize(N+1,0);
        dk[i].resize(N+1,0);
        lk[i].resize(N+1,0);
    }

    int u,v;
    for (int i=1;i<=NE;i++){
        cin>>u>>v;
        u++,v++;
        e[u].push_back(v);
        e[v].push_back(u);
        k[u]++,k[v]++;
        A[u][v]=1,A[v][u]=1;
    }

    heuristic();
    ACO();

    printf("\nTime taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

}
