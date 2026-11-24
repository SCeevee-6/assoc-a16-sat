// Consume the repository's sorted preorder codes and edge arrays. A triangle
// (a,k,b) corresponds to the internal node spanning leaves [a,b).
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
using Tri = array<uint16_t,14>;
using Adj = array<int,13>;
int ids[16][16][16], n, bits;
uint32_t code;
int pos, used;
int decode(int a, Tri& t) {
    if (pos > bits) throw runtime_error("bad tree");
    int p=pos++;
    if (p==bits || !(code & (1u<<(bits-1-p)))) return a+1;
    int k=decode(a,t), b=decode(k,t);
    if(b>=n || used>=n-2) throw runtime_error("bad span");
    t[used++]=ids[a][k][b];
    return b;
}
template<class T> void readn(ifstream& f,T* p,size_t count){
    if(!f.read(reinterpret_cast<char*>(p),sizeof(T)*count)) throw runtime_error("short graph");
}
int main(int argc,char**argv){
    if(argc!=2)return 2;
    ifstream f(argv[1],ios::binary);
    int32_t nv; int64_t ne; readn(f,&nv,1); readn(f,&ne,1);
    n=nv==742900?15:16; bits=2*(n-2);
    if((n==16 && nv!=2674440)||ne!=int64_t(nv)*(n-3)/2) throw runtime_error("bad counts");
    int nt=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++)ids[i][j][k]=nt++;
    vector<Tri> ts(nv); uint32_t prev=0;
    for(auto &t:ts){readn(f,&code,1); if(code<=prev)throw runtime_error("duplicate or unordered code");prev=code;
        pos=used=0;if(decode(0,t)!=n-1||pos!=bits+1||used!=n-2)throw runtime_error("invalid tree");
        sort(t.begin(),t.begin()+used);}
    vector<int32_t> u(ne);readn(f,u.data(),ne);
    vector<Adj> adj(nv);vector<int> degree(nv);
    for(int64_t i=0;i<ne;i++){int v;readn(f,&v,1);int a=u[i];
        if(a<0||a>=v||v>=nv||degree[a]>=n-3||degree[v]>=n-3)throw runtime_error("bad edge");
        adj[a][degree[a]++]=v;adj[v][degree[v]++]=a;}
    for(int i=0;i<nv;i++){auto &row=adj[i];sort(row.begin(),row.begin()+degree[i]);
        if(degree[i]!=n-3||adjacent_find(row.begin(),row.begin()+degree[i])!=row.begin()+degree[i])throw runtime_error("bad neighbors");}
    string seed;vector<unsigned char> first(nv);vector<int> parent(nv),queue(nv);vector<int8_t> side(nv);
    while(cin>>seed){
        if(int(seed.size())!=nt||seed.find_first_not_of("01")!=string::npos)return 2;
        for(int i=0;i<nv;i++){int b=0;for(int j=0;j<n-2;j++)b^=seed[ts[i][j]]-'0';first[i]=b;}
        fill(side.begin(),side.end(),-1); int x=-1,y=-1;
        for(int root=0;root<nv&&x<0;root++)if(side[root]<0){
            int head=0,tail=0;queue[tail++]=root;side[root]=0;parent[root]=-1;
            while(head<tail&&x<0){int a=queue[head++];
                for(int j=0;j<degree[a];j++){int b=adj[a][j];if(first[a]!=first[b])continue;
                    if(side[b]<0){side[b]=1-side[a];parent[b]=a;queue[tail++]=b;}
                    else if(side[a]==side[b]){x=a;y=b;break;}}}}
        if(x<0){ofstream out("candidate.colors",ios::binary);for(int i=0;i<nv;i++)out.put(2*first[i]+side[i]);out.close();if(!out)return 3;cout<<"FOUND"<<endl;continue;}
        vector<int>a,b;for(int i=x;i>=0;i=parent[i])a.push_back(i);for(int i=y;i>=0;i=parent[i])b.push_back(i);
        int common=-1;while(!a.empty()&&!b.empty()&&a.back()==b.back()){common=a.back();a.pop_back();b.pop_back();}
        a.push_back(common);a.insert(a.end(),b.rbegin(),b.rend());
        if(a.size()%2==0)throw runtime_error("even cycle");
        cout<<"CYCLE";for(size_t i=0;i<a.size();i++){auto &s=ts[a[i]],&t=ts[a[(i+1)%a.size()]];array<int,28>d;
            auto end=set_symmetric_difference(s.begin(),s.begin()+n-2,t.begin(),t.begin()+n-2,d.begin());
            if(end-d.begin()!=4)throw runtime_error("nonflip edge");for(int j=0;j<4;j++)cout<<" "<<d[j];}
        cout<<endl;
    }
}
