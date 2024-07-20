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
#include <chrono>
#include <map>
#include <chrono>
using namespace std;
const int MAXN = 207, MAXM = 1000;
const int Q_LIM = 2007, K=40;
int SEARCH_WAY = 1,RANGE_LIM = 1, GRAPH = 1,VALUE_LIM = 20, REPEAT_PATH = 5;
//0 default 1 count_top 2 count_bottom

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
// int n=10,m=30;
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
        if (SEARCH_WAY % 4 == 1 or SEARCH_WAY % 4 == 3 and rand() % 2 == 0) {
            // sort(q + op, q + ed, cmp_1);
            for(int x = op; x < ed; x++) {
                if (q[x].second.vis.count() > q[op].second.vis.count()) swap(q[op],q[x]);
            }
        }
        if (SEARCH_WAY % 4 == 2 or SEARCH_WAY % 4 == 3 and rand() % 2 == 1) {
            // sort(q + op, q + ed, cmp_2);
             for(int x = op; x < ed; x++) {
                if (q[x].second.vis.count() < q[op].second.vis.count()) swap(q[op],q[x]);
            }
        }
        const auto &tp = q[op];
        int x = tp.first;
        // cout << "op" << op <<" " << "ed" <<" " << ed << " " << x <<" "<<tp.second.first<<endl;
        //顺序bug，艹
        map<int,int> histoty;
        for(auto it = to[x].begin();it != to[x].end() && ed < Q_LIM;it++) {
            int y = it->first;
            int pos = it->second;
            auto new_bit = (tp.second.channel) | h[pos];
            // cout <<"st:" <<x<<"tar:"<<y<<endl;
            // cout <<"vised:" <<tp.second.second<<"\ntp.second:"<<tp.second.first<<"h[pos]:"<<h[pos]<<endl;
            if(tp.second.vis[y]) continue;
            if(new_bit.all()) continue;
            if(histoty.count(y) != 0 and histoty[y] <= h[pos].count()) {
                continue;
            }
            else {
                histoty[y] = h[pos].count();
            }
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
    for(int i = 0;i < RANGE_LIM;i++) x[i] = rand() % (r - l + 1);
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
    tp.v = VALUE_LIM - (int)sqrt(rand() % (VALUE_LIM * VALUE_LIM));
    return;
}

void cover(const tag &tp) {
    bitset<K>st;
    for(int pid = tp.l;pid <= tp.r;pid++) {
        st.set(pid);
    }
    for(auto id : tp.path) h[id] |= st;
    return;
}

struct Point {
    double x, y;
};

double distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p2.y - p1.y, 2));
}

bool circle_contains_no_other_points(const vector<Point>& points, int i, int j) {
    const Point& u = points[i];
    const Point& v = points[j];
    Point midpoint = {(u.x + v.x) / 2.0, (u.y + v.y) / 2.0};
    double radius = distance(u, v) / 2.0;
    for (int k = 0; k < points.size(); ++k) {
        if (k != i && k != j) {
            if (distance(points[k], midpoint) <= radius) {
                return false;
            }
        }
    }
    return true;
}

vector<Point> generate_random_points(int n) {
    vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        points[i] = {static_cast<double>(rand()) / RAND_MAX, static_cast<double>(rand()) / RAND_MAX};
    }
    return points;
}

vector<Point> generate_points_boundary(int n, double margin = 0.1) {
    vector<Point> points(n);
    int idx = 0;
    for (int i = 0; i < n / 4; ++i) {
        points[idx++] = {static_cast<double>(rand()) / RAND_MAX, margin * static_cast<double>(rand()) / RAND_MAX};
        points[idx++] = {static_cast<double>(rand()) / RAND_MAX, 1.0 - margin * static_cast<double>(rand()) / RAND_MAX};
        points[idx++] = {margin * static_cast<double>(rand()) / RAND_MAX, static_cast<double>(rand()) / RAND_MAX};
        points[idx++] = {1.0 - margin * static_cast<double>(rand()) / RAND_MAX, static_cast<double>(rand()) / RAND_MAX};
    }
    return points;
}

vector<Point> generate_points_ring(int n, double radius = 0.5, double width = 0.1) {
    vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        double angle = 2 * M_PI * static_cast<double>(rand()) / RAND_MAX;
        double r = radius + (static_cast<double>(rand()) / RAND_MAX - 0.5) * width;
        points[i] = {r * cos(angle) + 0.5, r * sin(angle) + 0.5};
    }
    return points;
}

vector<Point> generate_points_clusters(int n, int num_clusters = 4, double cluster_std = 0.05) {
    vector<Point> points(n);
    vector<Point> centers(num_clusters);
    for (int i = 0; i < num_clusters; ++i) {
        centers[i] = {static_cast<double>(rand()) / RAND_MAX, static_cast<double>(rand()) / RAND_MAX};
    }
    for (int i = 0; i < n; ++i) {
        Point& center = centers[i % num_clusters];
        points[i] = {center.x + cluster_std * ((static_cast<double>(rand()) / RAND_MAX) * 2 - 1),
                     center.y + cluster_std * ((static_cast<double>(rand()) / RAND_MAX) * 2 - 1)};
    }
    return points;
}

vector<pair<int, int>> generate_gabriel_graph(const vector<Point>& points) {
    vector<pair<int, int>> edges;
    int n = points.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (circle_contains_no_other_points(points, i, j)) {
                edges.push_back({i + 1, j + 1});
            }
        }
    }
    return edges;
}


double get_cover(){
    double cnt = 0, sum = 0;
    for(int i = 0;i < m;i++) {
        sum += K;
        cnt += h[i].count();
    }
    return cnt / sum;
}


auto get_paths() {
    int g[MAXN][MAXN];
    memset(g,0x3f3f3f3f,sizeof(g));
    for(auto [x,y] : edge) {
        g[x][y] = g[y][x] = 1;
    }
    for(int k = 1;k <= n;k++) 
        for(int i = 1;i <= n;i++) 
            for(int j = 1;j <= n;j++) 
                g[i][j] = min(g[i][j],g[i][k] + g[k][j]);

    vector<pair<int,pair<int,int> >> res;
    for(int i = 1;i <= n;i++) {
        for(int j = i + 1;j <= n;j++) {
            res.push_back(make_pair(-g[i][j],make_pair(i,j)));
        }
    }
    sort(res.begin(),res.end());

    // res.erase(res.begin() + int(sqrt(res.size())),res.end());
    // random_shuffle(res.begin(),res.end());
    return res;
}

int get_longest() {
    int g[MAXN][MAXN];
    memset(g,0x3f,sizeof(g));
    for(auto [x,y] : edge) {
        g[x][y] = g[y][x] = 1;
    }
    for(int k = 1;k <= n;k++) 
        for(int i = 1;i <= n;i++) 
            for(int j = 1;j <= n;j++) 
                g[i][j] = min(g[i][j],g[i][k] + g[k][j]);

    int res = 0;
    for(int i = 1;i <= n;i++) 
        for(int j = 1;j <= n;j++)
            res = max(res,g[i][j]);
    return res;
}



int main(){
    int total_time = 0;
    int goon = true;
    while(goon) {
        
        cin >> n >> m;
        for(int i = 1;i <= n;i++) cin>>p[i];
        for(int i = 0,x,y;i < m;i++) {
            cin>>x>>y;
            edge.push_back(make_pair(x,y));
        }
        int J;
        cin >> J;
        cout <<"ASAS"<<endl;
        for(int i = 0;i < J;i++) {
            tag tp;
            int lim;
            cin >> tp.src >> tp.snk >> lim >> tp.l >> tp.r >> tp.v;
            tp.l--;
            tp.r--;
            while(lim--) {
                int y;
                cin >> y;
                tp.path.push_back(y);
                for(int pp = tp.l;pp <= tp.r;pp++) h[y-1].set(pp);
            }
            mission.push_back(tp);
        }

        cout <<"ASAS"<<endl;
        // for(int i = 0;i < edge.size();i++) {
        //     cout << i <<" " <<h[i]<<endl;
        // }
        int t = 50;
        cin >> t;
        cout <<"T: "<<t<<endl;
        int bt = t;
        int v = 0;
        while(t--) {
            while(1) {
                int x;
                cin >> x;
                v++;
                if (x == -1) {v--;break;}
            }
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
        goon = 0;cout<<"mission = "<<mission.size()<<endl;
        cout<<"cover ratio:"<<get_cover()<<endl;
        cout<<"longest"<<get_longest()<<endl;
        cout<<"ave j "<<double(v) / bt<<" " <<double(v) <<" "<< bt<<endl;
    }
    // cout << "1 case time" << total_time << endl;
    

    // for(auto [x,y] : edge) {
    //     cout << x <<" " << y << endl;
    // }
    return 0;
}

/*
RANGE_LIM 1 SEARCH_WAY 4 graph 0

*/