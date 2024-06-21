#include <bits/stdc++.h>
constexpr int k = 40;
constexpr int MAXN = 256;
constexpr int MAXM = 1024;
const int MAXJ = 6000;
int n, m, q, p[MAXN];
struct edge_t {
    uint64_t channel;
    std::vector<int> occupied;
    int first;
    int second;
    int index;
    bool deleted;
    void set(int l, int r) {
        uint64_t mask = (1ull << (r + 1)) - (1ull << l);
        #ifdef __SMZ_RUNTIME_CHECK
        if (l <= 0 || r > k || l > r) {
            abort();
        }
        #endif
        channel |= mask;
    }
    void clear(int l, int r) {
        uint64_t mask = (1ull << (r + 1)) - (1ull << l);
        #ifdef __SMZ_RUNTIME_CHECK
        if (l <= 0 || r > k || l > r) {
            abort();
        }
        #endif
        channel &= ~mask;
    }
    bool empty(int l, int r) {
        uint64_t mask = (1ull << (r + 1)) - (1ull << l);
        #ifdef __SMZ_RUNTIME_CHECK
        if (l <= 0 || r > k || l > r) {
            abort();
        }
        #endif
        return (channel & mask) == 0;
    }
    void remove(int x) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (occupied.empty() || x <= 0 || x > q) {
            abort();
        }
        #endif
        for (int i = 0; i < occupied.size(); ++i) {
            if (occupied[i] == x) {
                std::swap(occupied[i], occupied.back());
                occupied.pop_back();
                return;
            }
        }
        #ifdef __SMZ_RUNTIME_CHECK
        abort();
        #endif
    }
    void insert(int x) {
        #ifdef __SMZ_RUNTIME_CHECK
        int cnt = 0;
        for (int i = 0; i < occupied.size(); ++i) {
            if (occupied[i] == x) {
                cnt += 1;
            }
        }
        if (cnt > 1 || x <= 0 || x > q) {
            abort();
        }
        #endif
        occupied.push_back(x);
    }
} edges[MAXM];
struct query_t {
    int from;
    int to;
    int l;
    int r;
    int value;
    int index;
    int span;
    bool dead;
    std::vector<std::pair<int, int>> path, backup;
    void apply(const std::vector<std::pair<int, int>>& new_path) const {
        int node = from;
        int channel = -1;
        for (auto [e, L] : new_path) {
            edges[e].set(L, L + span);
            edges[e].insert(index);
            if (channel != -1 && channel != L) {
                p[node] -= 1;
            }
            #ifdef __SMZ_RUNTIME_CHECK
            if (node != edges[e].first && node != edges[e].second) {
                abort();
            }
            if (p[node] < 0 || edges[e].deleted) {
                abort();
            }
            #endif
            node = (node != edges[e].first ? edges[e].first : edges[e].second);
            channel = L;
        }
    }
    void undo(const std::vector<std::pair<int, int>>& the_path) const {
        int node = from;
        int channel = -1;
        for (auto [e, L] : the_path) {
            edges[e].clear(L, L + span);
            edges[e].remove(index);
            if (channel != -1 && channel != L) {
                p[node] += 1;
            }
            #ifdef __SMZ_RUNTIME_CHECK
            if (node != edges[e].first && node != edges[e].second) {
                abort();
            }
            #endif
            node = (node != edges[e].first ? edges[e].first : edges[e].second);
            channel = L;
        }
    }
    void undo() {
        backup = std::move(path);
        path.clear();
        undo(backup);
    }
    void redo() {
        path = backup;
        apply(path);
    }
    void plan(const std::vector<std::pair<int, int>>& new_path) {
        path = new_path;
        apply(path);
        apply(backup);
    }
    void confirm() {
        undo(backup);
        undo(path);
        apply(path);
    }
} query[MAXJ];
std::vector<std::pair<int, edge_t*>> G[MAXN];
namespace io {
    constexpr int MAXBUFFER = 1024 * 1024 * 8;
    char ibuffer[MAXBUFFER], *iptr, obuffer[MAXBUFFER], *optr;
    inline void start_reading() { // 开始读取新的一行
        std::ignore = fgets(ibuffer, sizeof(ibuffer), stdin);
        iptr = ibuffer;
    }
    inline void start_writing() { //开始输出新的一行
        optr = obuffer;
    }
    inline int read_int() { //读入有符号整数
        char* nxt;
        int ret = strtol(iptr, &nxt, 10);
        iptr = nxt;
        return ret;
    }
    inline double read_double() noexcept { // 读入浮点数
        char* nxt;
        double ret = strtod(iptr, &nxt);
        iptr = nxt;
        return ret;
    }
    inline void write_int(int val) { //输出有符号整数，输出完一行后需要调用flush。
        char tmp[32], *now = tmp + 20;
        int length = 1;
        if (val < 0) {
            *optr++ = '-';
            val *= -1;
        }
        *now = ' ';
        do {
            *--now = '0' + val % 10;
            val /= 10;
            length += 1;
        } while (val > 0);
        memcpy(optr, now, length);
        optr += length;
    }
    inline void newline() {
        optr[-1] = '\n';
    }
    inline void flush() {
        if (optr != obuffer) {
            optr[-1] = '\n';
        }
        fwrite(obuffer, 1, optr - obuffer, stdout);
        fflush(stdout);
    }
}
namespace testcase {
    int n, m, q, p[MAXN];
    std::pair<int, int> nodes[MAXM];
    std::vector<std::tuple<int, int, int, int, int, int, std::vector<int>>> business;
    void run() {
        io::start_reading();
        n = io::read_int();
        m = io::read_int();
        io::start_reading();
        for (int i = 1; i <= n; ++i) {
            p[i] = io::read_int();
        }
        for (int i = 1; i <= m; ++i) {
            io::start_reading();
            int x = io::read_int();
            int y = io::read_int();
            nodes[i] = {x, y};
        }
        io::start_reading();
        q = io::read_int();
        for (int j = 1; j <= q; ++j) {
            io::start_reading();
            int s = io::read_int();
            int t = io::read_int();
            int m = io::read_int();
            int L = io::read_int();
            int R = io::read_int();
            int V = io::read_int();
            io::start_reading();
            std::vector<int> E;
            for (int i = 1; i <= m; ++i) {
                int e = io::read_int();
                E.push_back(e);
            }
            business.emplace_back(s, t, m, L, R, V, std::move(E));
        }
    }
    void start() {
        ::n = n;
        ::m = m;
        for (int i = 1; i <= n; ++i) {
            ::p[i] = p[i];
        }
        for (int i = 1; i <= n; ++i) {
            ::G[i].clear();
        }
        for (int i = 1; i <= m; ++i) {
            auto [x, y] = nodes[i];
            ::edges[i] = {};
            ::edges[i].first = x;
            ::edges[i].second = y;
            ::edges[i].index = i;
            ::edges[i].deleted = false;
            ::G[x].emplace_back(y, &::edges[i]);
            ::G[y].emplace_back(x, &::edges[i]);            
        }
        ::q = q;
        for (int i = 1; i <= business.size(); ++i) {
            const auto& [s, t, _, L, R, V, E] = business[i - 1];
            ::query[i].from = s;
            ::query[i].to = t;
            ::query[i].l = L;
            ::query[i].r = R;
            ::query[i].value = V;
            ::query[i].index = i;
            ::query[i].span = R - L;
            ::query[i].dead = false;
            ::query[i].path.clear();
            for (auto e : E) {
                ::query[i].path.emplace_back(e, L);
            }
            ::query[i].apply(::query[i].path);
        }
    }
}
std::vector<std::pair<int, int>> bfs(const query_t& qry) {
    static int vis[MAXN];
    static std::pair<int, int> father[MAXN];
    for (int i = 1; i <= n; ++i) {
        vis[i] = false;
    }
    std::queue<std::pair<int, int>> queue;
    queue.emplace(qry.from, qry.l);
    while (!queue.empty()) {
        auto [x, i] = queue.front();
        queue.pop();
        if (x == qry.to) {
            break;
        }
        for (auto [y, info] : G[x]) {
            if (!info->empty(i, i + qry.span)) {
                continue;
            }
            if (!vis[y]) {
                vis[y] = true;
                queue.emplace(y, i);
                father[y] = {x, info->index};
            }
        }
    }
    if (!vis[qry.to]) {
        return {};
    }
    std::vector<std::pair<int, int>> path;
    int node = qry.to;
    while (node != qry.from) {
        auto [prev, e] = father[node];
        path.emplace_back(e, qry.l);
        node = prev;
    }
    std::reverse(path.begin(), path.end());
    return path;
}
std::vector<int> solve(int e) {
    int s = edges[e].first, t = edges[e].second;
    for (int i = 0; i < G[s].size(); ++i) {
        if (G[s][i].second->index == e) {
            G[s].erase(G[s].begin() + i);
            break;
        }
    }
    std::swap(s, t);
    for (int i = 0; i < G[s].size(); ++i) {
        if (G[s][i].second->index == e) {
            G[s].erase(G[s].begin() + i);
            break;
        }
    }
    std::vector<int> ret;
    std::vector<int> deleted = edges[e].occupied;
    #ifdef __SMZ_RUNTIME_CHECK
    for (int i = 0; i < deleted.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            if (deleted[i] == deleted[j]) {
                abort();
            }
        }
    }
    #endif
    for (auto i : deleted) if (!query[i].dead) {
        query[i].undo();
        auto new_path = bfs(query[i]);
        if (new_path.empty()) {
            query[i].redo();
            query[i].dead = true;
        }
        else {
            query[i].plan(new_path);
            ret.push_back(i);
        }
    }
    for (auto i : ret) {
        query[i].confirm();
    }
    return ret;
}
int main() {
    #ifdef __SMZ_RUNTIME_CHECK
    std::ignore = freopen("../release/testcase1.in", "r", stdin);
    std::ignore = freopen("../release/output.txt", "w", stdout);
    #endif
    testcase::run();
    io::start_reading();
    int T = io::read_int();
    while (T--) {
        testcase::start();
        for (;;) {
            io::start_reading();
            int e = io::read_int();
            if (e == -1) {
                break;
            }
            auto indices = solve(e);
            io::start_writing();
            io::write_int(indices.size());
            io::newline();
            for (auto i : indices) {
                io::write_int(i);
                io::write_int(query[i].path.size());
                io::newline();
                for (auto [e, c] : query[i].path) {
                    io::write_int(e);
                    io::write_int(c);
                    io::write_int(c + query[i].span);
                }
                io::newline();
            }
            io::flush();
            edges[e].deleted = true;
        }
    }
    return 0;
}
