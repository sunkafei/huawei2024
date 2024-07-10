#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <bitset>
#include <queue>
#include <cstring>
#include <ctime>
#include <fstream>
#include <ostream>
#include <set>
#include <assert.h>
#include <cmath>
using namespace std;
const int MAXN = 207, MAXM = 1000;
const int Q_LIM = 2007, K=40;



struct tag{
    int src,snk;
    int s,l,r,v;
    vector<int>path;
};
vector<pair<int,int> >edge;
vector<pair<int,int> >to[MAXN];
bitset<K> h[MAXM];
vector<tag>mission;
int n = 200,m = 1000;
int p[207] = {0};



int op = 0,ed = 1;
pair<int,pair<bitset<K>,bitset<MAXN> > > q[Q_LIM];
int from[Q_LIM] = {0},from_edge[Q_LIM] = {0};


pair<int,int> findZeroSegments(const bitset<K>& bs) {
    vector<pair<int,int>> segments;
    bool inSegment = false;
    int start = 0;

    for (int i = 0; i < bs.size(); ++i) {
        if (!bs[i] && !inSegment) {
            inSegment = true;
            start = i;
        } else if (bs[i] && inSegment) {
            segments.push_back({start, static_cast<int>(i - 1)});
            inSegment = false;
        }
    }

    if (inSegment) {
        segments.push_back({start, static_cast<int>(bs.size() - 1)});
    }

    int x = rand() % segments.size();
    return segments[x];
}

bool cmp_1(const pair<int,pair<bitset<K>,bitset<MAXN> > > &a,const pair<int,pair<bitset<K>,bitset<MAXN> > > &b) {
    return a.second.second.count() > b.second.second.count();
}

bool cmp_2(const pair<int,pair<bitset<K>,bitset<MAXN> > > &a,const pair<int,pair<bitset<K>,bitset<MAXN> > > &b) {
    return a.second.second.count() < b.second.second.count();
}

bool get_path(int s,int t) {
    bitset<K> blank_1;
    bitset<MAXN> blank_2;
    blank_2.set(s);
    op = 0,ed = 1;
    q[0] = make_pair(s,make_pair(blank_1,blank_2));
    memset(from,0,sizeof(from));
    // cout <<"--------------------------------"<<endl;
    // cout << "get_path" <<s <<" "<<t<<endl;

    while(op < ed) {
        int r = rand();
        if (r % 2 == 0) sort(q + op, q + ed, cmp_1);
        else sort(q + op, q + ed, cmp_2);
        auto tp = q[op];
        int x = tp.first;
        // cout << "op" << op <<" " << "ed" <<" " << ed << " " << x <<" "<<tp.second.first<<endl;
        //顺序bug，艹
        for(auto it = to[x].begin();it != to[x].end() && ed < Q_LIM;it++) {
            int y = it->first;
            int pos = it->second;
            auto new_bit = (tp.second.first) | h[pos];
            // cout <<"st:" <<x<<"tar:"<<y<<endl;
            // cout <<"vised:" <<tp.second.second<<"\ntp.second:"<<tp.second.first<<"h[pos]:"<<h[pos]<<endl;
            if(tp.second.second[y]) continue;
            if(new_bit.all()) continue;
            bitset<MAXN> new_vis = tp.second.second;
            new_vis.set(y);
            q[ed] = make_pair(y,make_pair(new_bit,new_vis));
            from[ed] = op;
            from_edge[ed] = pos;
            ed++;
            if(y == t) {
                op = ed;
                return true;
            }
        }
        op++;
    }
    return false;
}

pair<int,int> get_rand_seg(pair<int,int> seg){
    int l = seg.first;
    int r = seg.second;
    int x[52] = {0};
    for(int i = 0;i < 50;i++) x[i] = rand() % (r - l + 1);
    // return make_pair(min(x,min(y,z)),max(x,(max(y,z))));
    // return make_pair(min(x,y),max(x,y));
    return make_pair(l , l + *max_element(x,x + 50));
}


void get_output(tag &tp) {
    tp.path.clear();
    // cout << "Q[op]:"<<q[op-1].first<<" "<<q[op-1].second.first << endl;
    pair<int,int> seg = findZeroSegments(q[op-1].second.first);
    seg = get_rand_seg(seg);
    bitset<K>st;
    for(int i = seg.first;i <= seg.second;i++) {
        st.set(i);
    }
    // cout <<"Q:"<<endl;
    // cout << "st:" << st << endl;
    // for(int i = 0;i < ed;i++) cout << q[i].first <<" ";cout << endl;
    // for(int i = 0;i < ed;i++) cout << from[i] <<" ";cout << endl;
    int pos = op-1;
    while(pos != 0) {
        int fr = q[from[pos]].first;
        int id = from_edge[pos];
        // cout <<"before" << endl;
        // cout << pos <<" "<< id <<" "<<h[id] << " "<< st<<endl;
        h[id] |= st;
        // cout <<"after" << endl;
        // cout << pos <<" "<<id <<" "<<h[id] << " "<< st<<endl;
        tp.path.push_back(id);
        pos = from[pos];
    }
    reverse(tp.path.begin(),tp.path.end());
    // cout << "from" << tp.src <<" to" << tp.snk << endl<<"id:";
    // for(auto x :tp.path){
    //     cout <<x << " ";   
    // }cout <<endl;
    tp.s = tp.path.size();
    tp.l = seg.first;
    tp.r = seg.second;
    tp.v = rand() % 10;
    return;
}


bool check(tag x) {
    int s = x.src;
    for(auto i : x.path) {
        int a,b;
        a = edge[i].first;
        b = edge[i].second;
        if(s == a) s = b;
        else if (s == b)
        {
            s = a;
        }
        else {
            cout<<"gg";
            return false;
        }
        
    }
    return true;
}
int main(){
    cin>>n>>m;
    for(int i = 1;i <= n;i++) cin>>p[i];
    edge.push_back({0,0});
    for(int x,y,i = 0;i < m;i++) {
        cin >> x >> y;
        edge.push_back(make_pair(x,y));
    }
    int J;
    cin >> J;
    for(int l,i = 0;i < J;i++) {
        tag x;
        cin >> x.src >> x.snk >> l >> x.l >> x.r >> x.v;
        for(int j = 0;j < l;j++) {
            int sd;
            cin >> sd;
            x.path.push_back(sd);
        }
        if(!check(x)) {
            cout << x.path.size() << endl;
            cout << x.src <<" " << x.snk << endl;
            for(auto gg : x.path) {
                cout << edge[gg].first << " " << edge[gg].second << endl;
            }
            cout <<" -------------------- " << endl;
        }
        mission.push_back(x);
    }
    // for(int i = 0;i < edge.size();i++) {
    //     cout << i <<" " <<h[i]<<endl;
    // }
    // t = 40;
    // while(t) {
    //     int x = 1;
    //     int kill[MAXM] = {0};
    //     for(int j = 1;j <= m;j++) kill[j] = j;
    //     random_shuffle(kill+1,kill + m + 1);
    //     for(int i = 1;i <= min(x,m);i++) outfile << kill[i] << endl;
    //     outfile << -1 << endl;
    //     t--;
    // }
    return 0;
}