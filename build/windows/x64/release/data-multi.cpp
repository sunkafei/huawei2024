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
using namespace std;
const int MAXN = 207, MAXM = 1000;
const int Q_LIM = 2007, K=40;
int SEARCH_WAY = 1,RANGE_LIM = 2, GRAPH = 1,VALUE_LIM = 30, REPEAT_PATH = 20;
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
int gp[207][207] = {0};


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
        // random_shuffle(to[x].begin(),to[x].end());
        sort(to[x].begin(),to[x].end(),[t](auto a, auto b) {
        return gp[a.first][t] > gp[b.first][t];});
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
    int delta = 0;
    while(rand() % 2 == 0) {
        delta++;
    }
    delta = min(delta,r - l);
    // return make_pair(min(x,min(y,z)),max(x,(max(y,z))));
    // return make_pair(min(x,y),max(x,y));
    return make_pair(l , l + delta);
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
        for(int j = 1;j <= n;j++) {
            res = max(res,g[i][j]);
            gp[i][j] = g[i][j];
        }
    return res;
}



int main(){
    int total_time = 0;
    int goon = true;
    while(goon) {
        srand(time(NULL));
        ofstream outfile;
        // string rule;
        // cin >> rule;
        outfile.open("smz.in");

        for(int i = 0;i <= n;i++) {
            p[i] = rand() % 16 + 5;
            if (rand() % 7 == 0) p[i] = 0;
            to[i].clear();
        }
        edge.clear();
        mission.clear();
        for(int i = 0;i <= m;i++) h[i].reset();


        if (!GRAPH) {//不生成gg图
                for(int i = 1;i < n;i++) {
                edge.push_back(make_pair(i,i + 1));
            }
        }
        else {//尝试生成gg图
            vector<Point> points;
                switch (GRAPH) {
                case 1:
                    points = generate_random_points(n);
                    break;
                case 2:
                    points = generate_points_boundary(n);
                    break;
                case 3:
                    points = generate_points_ring(n);
                    break;
                case 4:
                    points = generate_points_clusters(n);
                    break;
                default:
                    cout << "Invalid ID!" << endl;
            }
            edge = generate_gabriel_graph(points);
        }
        cout << "init edge.size() " << edge.size() << endl;
        // m = edge.size();
        while(edge.size() < m) {
            int x = rand() % n + 1;
            int y = rand() % n + 1;
            while(y == x) y = rand() % n + 1;
            if(x > y) swap(x,y);
            // x = 1;
            // y = n;   
            if (find(edge.begin(),edge.end(),make_pair(x,y)) == edge.end()) {
                continue;
            }
            edge.push_back(make_pair(x,y));
        }
        cout << "total edge.size() " << edge.size() << endl;
        random_shuffle(edge.begin(),edge.end());
        for(int i = 0;i < edge.size();i++) {
            to[edge[i].first].push_back(make_pair(edge[i].second,i));
            to[edge[i].second].push_back(make_pair(edge[i].first,i));
        }
        cout<<"longest"<<get_longest()<<endl;
        // set<pair<int,int>> S;
        // for(auto x : edge) S.insert(make_pair(min(x.first,x.second),max(x.first,x.second)));
        // cout << S.size() <<" " << edge.size()<<endl;
        // assert(S.size() == edge.size());
        // sort(edge.begin(),edge.end());


        auto paths = get_paths();
        // paths.erase(paths.begin() + (int)sqrt(paths.size()),paths.end());
        paths.erase(paths.begin() + paths.size() / 50,paths.end());
        random_shuffle(paths.begin(),paths.end());
        cout <<"get paths done "<< endl;
        int minn = m,maxx = 0;
        for(auto [x,y] : paths) {
            static int pathi = 0;
            pathi++;
            if(pathi % 100 == 0) cout << "part:" << pathi << " " << "(mission)" <<" " << mission.size() <<"/" << paths.size()<< endl;
            tag tp;
            tp.src = y.first;
            tp.snk = y.second;
            if(get_path(tp.src,tp.snk)) {
                get_output(tp);
                int T = REPEAT_PATH;
                while(T--) {
                    tag new_tp;
                    new_tp.src = tp.src;
                    new_tp.snk = tp.snk;
                    if(get_path(new_tp.src,new_tp.snk)) {
                        get_output(new_tp);
                        if(new_tp.path.size() > tp.path.size()) tp = new_tp;
                    }
                }
                // cout <<" pre success" << tp.src <<" "<<tp.snk << " " << mission.size() <<endl;
                cover(tp);
                // if(tp.path.size() >= 30) continue;
                mission.push_back(tp);
                minn = min((int)tp.path.size(),minn);
                maxx = max((int)tp.path.size(),minn);
                // cout <<" las success" << tp.src <<" "<<tp.snk << " " <<mission.size() <<endl;
            }
        }
        cout<<"init mission = "<<mission.size()<<endl;
        cout<<"init cover ratio:"<<get_cover()<<endl;
        cout<<"init min dis:"<<minn<<endl;
        cout<<"init max dis:"<<maxx<<endl;
        int J = 5000;

        for(int i = 1;i <= 20000 and mission.size() < J;i++) {
            if(i % 1000 == 0) cout << "pre part:" << i << " " << "(mission)" <<" " << mission.size() << " " << mission.back().path.size() << endl;
            tag tp;
            tp.src = rand() % n + 1;
            tp.snk = rand() % n + 1;
            while(tp.src == tp.snk) tp.snk = rand() % n + 1;
            if(get_path(tp.src,tp.snk)) {
                get_output(tp);
                int T = REPEAT_PATH;
                while(T--) {
                    tag new_tp;
                    new_tp.src = tp.src;
                    new_tp.snk = tp.snk;
                    if(get_path(new_tp.src,new_tp.snk)) {
                        get_output(new_tp);
                        if(new_tp.path.size() > tp.path.size()) tp = new_tp;
                    }
                }
                if(tp.path.size() <= 4) continue;
                if(tp.path.size() <= 7 and rand() % 2 == 0) continue; 
                if(tp.path.size() <= 10 and rand() % 3 == 0) continue; 
                cover(tp);
                mission.push_back(tp);
                // cout <<" las success" << tp.src <<" "<<tp.snk << " " <<mission.size() <<endl;
            }
        }



        for(int i = 1;i <= 1000000 and mission.size() < J;i++) {
            if(i % 100000 == 0) cout << "lat part:" << i << " " << "(mission)" <<" " << mission.size() << endl;
            tag tp;
            tp.src = rand() % n + 1;
            tp.snk = rand() % n + 1;
            while(tp.src == tp.snk) tp.snk = rand() % n + 1;
            if(get_path(tp.src,tp.snk)) {
                get_output(tp);
                int T = REPEAT_PATH;
                while(T--) {
                    tag new_tp;
                    new_tp.src = tp.src;
                    new_tp.snk = tp.snk;
                    if(get_path(new_tp.src,new_tp.snk)) {
                        get_output(new_tp);
                        if(new_tp.path.size() > tp.path.size()) tp = new_tp;
                    }
                }
                // cout <<" pre success" << tp.src <<" "<<tp.snk << " " << mission.size() <<endl;
                // get_output(tp);
                // if(tp.path.size() <= 4) continue;
                // if(tp.path.size() <= 7 and rand() % 2 == 0) continue; 
                // if(tp.path.size() <= 10 and rand() % 3 == 0) continue; 
                cover(tp);
                mission.push_back(tp);
                // cout <<" las success" << tp.src <<" "<<tp.snk << " " <<mission.size() <<endl;
            }
        }
        cout<<"mission = "<<mission.size()<<endl;
        cout<<"cover ratio:"<<get_cover()<<endl;
        outfile<<n<<" "<<m<<endl;
        for(int i = 1;i <= n;i++) outfile<<p[i]<<" ";outfile<<endl;
        for(auto x : edge) {
            assert(x.first !=0 and x.second != 0);
            if (rand() & 1) outfile << x.first<<" " << x.second<<endl;
            else outfile << x.second<<" " << x.first<<endl;
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
        int t = 50;
        outfile << t << endl;
        
        double num[107] = {0};
        // for(int i = 1;i <= t;i++) num[i] = i;
        for(int i = 1;i <= t;i++) num[i] = rand();
        double sum = accumulate(num + 1,num + t + 1 ,0);

        int cnt[1007] = {0};
        for(auto x : mission) {
            for(auto y : x.path) {
                cnt[y]++;
            }
        }

        while(t) {
            int x = round(num[t] * 50 * t / sum);
            x = 60;
            int kill[MAXM] = {0};
            for(int j = 1;j <= m;j++) kill[j] = j;
            sort(kill + 1,kill + m + 1,[cnt](const int &a,const int &b){
                if(cnt[a] != cnt[b]) return cnt[a] > cnt[b];
                if(h[a].count() != h[b].count()) return h[a].count() > h[b].count();
                return a > b;
            });
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
        outfile.close();

        static int idx = 0;
        idx++;
        string command = "main\n";
        cout << "building data done. NO." << idx << endl;
        // auto start = std::chrono::high_resolution_clock::now();
        int res = system(command.c_str());
        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // int res = 0;
        goon = 0;
    }
    // cout << "1 case time" << total_time << endl;
    cout<<"mission = "<<mission.size()<<endl;
    cout<<"cover ratio:"<<get_cover()<<endl;
    cout<<"longest"<<get_longest()<<endl;

    // for(auto [x,y] : edge) {
    //     cout << x <<" " << y << endl;
    // }
    return 0;
}

/*
RANGE_LIM 1 SEARCH_WAY 4 graph 0

*/