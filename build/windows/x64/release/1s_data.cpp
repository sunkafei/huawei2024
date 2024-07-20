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
int SEARCH_WAY = 0,RANGE_LIM = 2;


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
// int n=30,m=100;
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
        if (SEARCH_WAY % 4 == 1 or SEARCH_WAY % 4 == 3 and rand() % 2 == 0) sort(q + op, q + ed, cmp_1);
        if (SEARCH_WAY % 4 == 2 or SEARCH_WAY % 4 == 3 and rand() % 2 == 1) sort(q + op, q + ed, cmp_2);
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
    tp.v = rand() % 50;
    return;
}


struct Point {
    double x, y;
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
    bool operator <(const Point& other) const {
        return x < other.x || x == other.x && y < other.y;
    }
};

double distanceSquared(const Point& p1, const Point& p2) {
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

bool areCollinear(const Point& p1, const Point& p2, const Point& p3) {
    return std::fabs((p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y)) < std::numeric_limits<double>::epsilon();
}

bool isGabrielEdge(const Point& p1, const Point& p2, const Point& p3) {
    double d12 = distanceSquared(p1, p2);
    double d13 = distanceSquared(p1, p3);
    double d23 = distanceSquared(p2, p3);
    // cout << "isGabrielEdge" << d12 << " " << d13 <<" "<< d23 <<" " << areCollinear(p1, p2, p3)<<endl;
    return d12 < (d13 + d23) && !areCollinear(p1, p2, p3);
}

std::vector<std::pair<int, int>> generateGabrielGraph(const std::vector<Point>& points) {
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < points.size(); ++i) {
        for (int j = i + 1; j < points.size(); ++j) {
            bool isGabriel = true;
            for (int k = 0; k < points.size(); ++k) {
                if (k != i && k != j && !isGabrielEdge(points[i], points[j], points[k])) {
                    isGabriel = false;
                    break;
                }
            }
            if (isGabriel) {
                edges.emplace_back(i + 1, j + 1);
            }
        }
    }
    return edges;
}

bool isDuplicate(const Point& a, const Point& b, double tolerance = 1e-6) {
    return std::fabs(a.x - b.x) < tolerance && std::fabs(a.y - b.y) < tolerance;
}


Point generateRandomPoint() {
    return {static_cast<double>(std::rand()) / RAND_MAX, 
            static_cast<double>(std::rand()) / RAND_MAX};
}

void GabrielGraph() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    std::set<Point> points;

    while (points.size() < n) {
        Point newPoint = generateRandomPoint();
        // std::cout << points.size() << endl;
        bool isDuplicateFound = false;
        for (const auto& existingPoint : points) {
            if (isDuplicate(newPoint, existingPoint)) {
                isDuplicateFound = true;
                break;
            }
        }
        // std::cout << "isDuplicateFound " << isDuplicateFound << endl;
        if (!isDuplicateFound) {
            bool isCollinearFound = false;
            for (const auto& otherPoint : points) {
                for (const auto& anotherPoint : points) {
                    if(otherPoint == anotherPoint) continue;
                    if (areCollinear(newPoint, otherPoint, anotherPoint)) {
                        isCollinearFound = true;
                        break;
                    }
                }
                if (isCollinearFound) break;
            }
            
            if (!isCollinearFound) {
                points.insert(newPoint);
            }
        }
    }
    // cout <<"points "<<points.size()<<endl;
    auto ggedges = generateGabrielGraph(std::vector<Point>(points.begin(),points.end()));
    edge.insert(edge.end(),ggedges.begin(),ggedges.end());    
    return;
}
int main(){
    int goon = true;
    while(goon) {
        srand(time(NULL));
        SEARCH_WAY = 0;
        RANGE_LIM = 1;
        ofstream outfile;
        // string rule;
        // cin >> rule;
        outfile.open("smz.in");

        for(int i = 0;i <= n;i++) {
            p[i] = rand() % 10 + 10;
            to[i].clear();
        }
        edge.clear();
        mission.clear();
        for(int i = 0;i <= m;i++) h[i].reset();


        int graph = 0;
        if (graph) {//不生成gg图
                for(int i = 1;i < n;i++) {
                edge.push_back(make_pair(i,i + 1));
            }
        }
        else {//尝试生成gg图
            GabrielGraph();
            // for(auto [x,y] : edge) {
            //     cout << x <<" " << y << endl;
            // }
        }
        cout <<"SEARCH_WAY: " << SEARCH_WAY << endl;
        cout <<"RANGE_LIM: " << RANGE_LIM << endl;
        cout <<"graph: " << graph << endl;

        cout << "init edge.size() " << edge.size() << endl;
        while(edge.size() < m) {
            int x = rand() % n + 1;
            int y = rand() % n + 1;
            while(y == x) y = rand() % n + 1;
            if(x > y) swap(x,y);
            if (find(edge.begin(),edge.end(),make_pair(x,y)) != edge.end()) {
                continue;
            }
            edge.push_back(make_pair(x,y));
        }
        cout << "total edge.size() " << edge.size() << endl;
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
        for(int i = 1;i <= 10000000 and mission.size() < J;i++) {
            // if(i % 10 == 0) cout << "edge:" << i << endl;
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
        cout<<"j="<<mission.size()<<endl;
        outfile<<n<<" "<<m<<endl;
        for(int i = 1;i <= n;i++) outfile<<p[i]<<" ";outfile<<endl;
        for(auto x : edge) {
            assert(x.first !=0 and x.second != 0);
            outfile << x.first<<" " << x.second<<endl;
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

        while(t) {
            int x = 60;
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
        outfile.close();

        static int idx = 0;
        idx++;
        string command = "main\n";
        cout << "building data done. NO." << idx << endl;
        int res = system(command.c_str());
        // int res = 0;
        goon = 0;
    }
    return 0;
}