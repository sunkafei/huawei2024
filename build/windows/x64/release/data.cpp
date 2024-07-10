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
#include <cstdlib>
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
struct q_element {
    bitset<K> channel;
    bitset<MAXN> vis;
    int from,from_edge;
    q_element(){};
    q_element(bitset<K> channel,bitset<MAXN> vis,int from,int from_edge):channel(channel),vis(vis),from(from),from_edge(from_edge){};
};
pair<int,q_element> q[Q_LIM];


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

bool cmp_1(const pair<int,q_element> &a,const pair<int,q_element> &b) {
    return a.second.vis.count() > b.second.vis.count();
}

bool cmp_2(const pair<int,q_element> &a,const pair<int,q_element> &b) {
    return a.second.vis.count() < b.second.vis.count();
}

bool get_path(int s,int t) {
    bitset<K> blank_1;
    bitset<MAXN> blank_2;
    blank_2.set(s);
    op = 0,ed = 1;
    q[0] = make_pair(s,q_element(blank_1,blank_2,-1,-1));
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
            auto new_bit = (tp.second.channel) | h[pos];
            // cout <<"st:" <<x<<"tar:"<<y<<endl;
            // cout <<"vised:" <<tp.second.second<<"\ntp.second:"<<tp.second.first<<"h[pos]:"<<h[pos]<<endl;
            if(tp.second.vis[y]) continue;
            if(new_bit.all()) continue;
            bitset<MAXN> new_vis = tp.second.vis;
            new_vis.set(y);
            q[ed] = make_pair(y,q_element(new_bit,new_vis,op,pos));
            // from[ed] = op;
            // from_edge[ed] = pos;
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
    pair<int,int> seg = findZeroSegments(q[op-1].second.channel);
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
        int fr = q[q[pos].second.from].first;
        int id = q[pos].second.from_edge;
        // cout <<"before" << endl;
        // cout << pos <<" "<< id <<" "<<h[id] << " "<< st<<endl;
        h[id] |= st;
        // cout <<"after" << endl;
        // cout << pos <<" "<<id <<" "<<h[id] << " "<< st<<endl;
        tp.path.push_back(id);
        pos = q[pos].second.from;
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
int main(){
    srand(time(NULL));
    ofstream outfile;
    string output_file;
    cin >> output_file;
    outfile.open(output_file);

    for(int i = 0;i <= n;i++) p[i] = rand() % 10 + 10;
    for(int i = 1;i < n;i++) {
        edge.push_back(make_pair(i,i + 1));
    }
    for(int i = n;i <= m;i++) {
        int x = rand() % n + 1;
        int y = rand() % n + 1;
        while(y == x) y = rand() % n + 1;
        if(x > y) swap(x,y);
        if (find(edge.begin(),edge.end(),make_pair(x,y)) != edge.end()) {
            i--;
            continue;
        }
        edge.push_back(make_pair(x,y));
    }
    random_shuffle(edge.begin(),edge.end());
    for(int i = 0;i < edge.size();i++) {
        if (rand() & 1) swap(edge[i].first,edge[i].second);
        to[edge[i].first].push_back(make_pair(edge[i].second,i));
        to[edge[i].second].push_back(make_pair(edge[i].first,i));
    }

    set<pair<int,int>> S;
    for(auto x : edge) S.insert(make_pair(min(x.first,x.second),max(x.first,x.second)));
    // cout << S.size() <<" " << edge.size()<<endl;
    assert(S.size() == edge.size());
    // sort(edge.begin(),edge.end());


    int J = 5000;
    for(int i = 1;i <= J;i++) {
        if(i % 500 == 0) cout << i << endl;
        tag tp;
        tp.src = rand() % n + 1;
        tp.snk = rand() % n + 1;
        while(tp.src == tp.snk) tp.snk = rand() % n + 1;
        if(get_path(tp.src,tp.snk)) {
            // cout <<" pre success" << tp.src <<" "<<tp.snk << " " << mission.size() <<endl;
            get_output(tp);
            mission.push_back(tp);
            // cout <<" las success" << tp.src <<" "<<tp.snk << " " <<mission.size() <<endl;
        }
    }
    cout<<mission.size()<<endl;
    outfile<<n<<" "<<m<<endl;
    for(int i = 1;i <= n;i++) outfile<<p[i]<<" ";outfile<<endl;
    for(auto x : edge) {
        outfile << x.first <<" " << x.second<<endl;
    }
    outfile<<mission.size()<<endl;
    for(auto x: mission) {
        outfile << x.src <<" " << x.snk<<" "<<x.path.size()<< " " <<x.l + 1<<" " <<x.r + 1<<" "<<x.v<<endl;
        for(auto j : x.path) {
            outfile << j + 1<<" ";
        }outfile << endl;
    }
    // for(int i = 0;i < edge.size();i++) {
    //     cout << i <<" " <<h[i]<<endl;
    // }
    int t = 100;
    outfile << 100 << endl;
    
    double num[107] = {0};
    // for(int i = 1;i <= t;i++) num[i] = i;
    for(int i = 1;i <= t;i++) num[i] = rand();
    double sum = accumulate(num + 1,num + t + 1 ,0);

    while(t) {
        int x = round(num[t] * 6000 / sum);
        x = min(x,50);
        int kill[MAXM] = {0};
        for(int j = 1;j <= m;j++) kill[j] = j;
        random_shuffle(kill+1,kill + m + 1);
        for(int i = 1;i <= min(x,m);i++) outfile << kill[i] << endl;
        outfile << -1 << endl;
        t--;
    }
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

    string command = "main <" + output_file;
    cout << command;
    system(command.c_str());
    return 0;
}