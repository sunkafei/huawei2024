#include <cstdint>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <chrono>
#include <random>
#include <queue>
#include <cstring>
#include <unordered_set>
#include <bitset>
#include <cmath>
#include <queue>
constexpr int INF = 1 << 30;
constexpr int k = 40;
constexpr int MAXK = 44;
constexpr int MAXN = 256;
constexpr int MAXM = 1024;
constexpr int MAXQ = 6000;
constexpr int MAXTIME = 89;
constexpr int MAXT1 = 30;
constexpr int MAXC = 60;
constexpr double MAXJACCARD = 0.5;
// Hyperparameters --------------------------------
constexpr int DEGREE = 30;
constexpr int MAXGENTIME = 50;
constexpr double EXPAND_RATIO = 0.2;
constexpr double TRANSFORM_EXP = 2.0 / 3.0;
// ------------------------------------------------
int64_t iterations = 0;
int num_operations = INF;
int n, m, q, p[MAXN];
using path_t = std::vector<std::pair<int, int>>;
using case_t = std::vector<std::pair<double, int>>;
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
    bool empty(uint64_t mask) {
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
    path_t path, backup;
    static inline int64_t vis[MAXQ];
    static inline int64_t timestamp = 1;
    template<bool first=true, bool update=false> void apply(const path_t& new_path) const {
        int node = from;
        int channel = -1;
        for (auto [e, L] : new_path) {
            edges[e].set(L, L + span);
            if constexpr (update) {
                edges[e].insert(index);
            }
            if (channel != -1 && channel != L) {
                if constexpr (first) {
                    p[node] -= 1;
                    vis[node] = timestamp;
                }
                else {
                    if (vis[node] != timestamp) {
                        p[node] -= 1;
                    }
                }
            }
#ifdef __SMZ_RUNTIME_CHECK
            if (node != edges[e].first && node != edges[e].second) {
                abort();
            }
            if (edges[e].deleted) { //p[node] < -1 || 
                abort();
            }
#endif
            node = (node != edges[e].first ? edges[e].first : edges[e].second);
            channel = L;
        }
    }
    template<bool first=true, bool update=false> void undo(const path_t& the_path) const {
        int node = from;
        int channel = -1;
        for (auto [e, L] : the_path) {
            edges[e].clear(L, L + span);
            if constexpr (update) {
                edges[e].remove(index);
            }
            if (channel != -1 && channel != L) {
                if constexpr (first) {
                    p[node] += 1;
                    vis[node] = timestamp;
                }
                else {
                    if (vis[node] != timestamp) {
                        p[node] += 1;
                    }
                }
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
        path = std::move(backup);
        backup.clear();
        apply(path);
    }
    void confirm(path_t&& new_path) {
        path = std::move(new_path);
        timestamp += 1;
        apply<true>(path);
        apply<false>(backup);
    }
    auto cancel() {
        timestamp += 1;
        undo<true>(backup);
        undo<false>(path);
        apply(backup);
        auto ret = std::move(path);
        path = std::move(backup);
        return ret;
    }
    void replace(path_t&& new_path) {
        undo<true, true>(path);
        apply<true, true>(new_path);
        path = std::move(new_path);
    }
    void init(const path_t& the_path) {
        apply<true, true>(the_path);
    }
} query[MAXQ];
struct transaction_t {
    int64_t loss;
    std::vector<std::pair<int, path_t>> result;
    std::vector<int> dead;
    void commit() {
        for (auto iter = result.begin(); iter != result.end(); ++iter) {
            auto [i, new_path] = std::move(*iter);
            query[i].replace(std::move(new_path));
        }
        for (auto i : dead) {
            query[i].dead = true;
        }
        #ifdef __SMZ_RUNTIME_CHECK
        for (const auto& [i, path] : result) {
            if (query[i].dead) {
                abort();
            }
        }
        for (int i = 1; i <= q; ++i) if (!query[i].dead) {
            std::unordered_set<int> S;
            for (auto [e, c] : query[i].path) {
                S.insert(e);
                if (std::count(edges[e].occupied.begin(), edges[e].occupied.end(), i) != 1) {
                    abort();
                }
            }
            for (int e = 1; e <= m; ++e) if (!S.count(e)) {
                if (std::count(edges[e].occupied.begin(), edges[e].occupied.end(), i) > 0) {
                    abort();
                }
            }
        }
        for (int i = 1; i <= n; ++i) if (p[i] < 0) {
            abort();
        }
        #endif
    }
};
struct item_t {
    double value;
    int counter;
    case_t sequence;
    item_t(double val, int cnt, case_t&& seq) : value(val), counter(cnt), sequence(std::move(seq)) {}
    bool operator< (const item_t& rhs) const {
        return value < rhs.value;
    }
};
std::mt19937 engine;
std::vector<std::pair<int, edge_t*>> G[MAXN];
std::vector<std::vector<std::pair<double, int>>> pretests;
const auto start_time = std::chrono::steady_clock::now();
template<typename... T> void print(const T&... sth) {
#ifdef __SMZ_NATIVE_TEST
    (..., (std::cerr << sth << " ")) << std::endl;
#endif
}
template<typename T> void print(const std::vector<T>& vec) {
    for (const auto& v : vec) {
        std::cerr << v << ' ';
    }
    std::cerr << std::endl;
}
template<typename T, int maxsize> class deque_t {
private:
    T* data;
    int L = maxsize;
    int R = maxsize;
public:
    deque_t() : data(new T[maxsize * 2]) {}
    T& front() {
        return data[L];
    }
    T& back() {
        return data[R - 1];
    }
    int size() const {
        return R - L;
    }
    bool empty() const {
        return L == R;
    }
    void clear() {
        L = maxsize;
        R = maxsize;
    }
    void push_front(const T& value) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (L <= 0) {
            abort();
        }
        #endif
        L -= 1;
        data[L] = value;
    }
    void push_back(const T& value) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R >= maxsize * 2) {
            abort();
        }
        #endif
        data[R] = value;
        R += 1;
    }
    void pop_front() {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R <= L) {
            abort();
        }
        #endif
        L += 1;
    }
    void pop_back() {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R <= L) {
            abort();
        }
        #endif
        R -= 1;
    }
    template<typename... F> void emplace_back(F&&... args) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R >= maxsize * 2) {
            abort();
        }
        #endif
        new(&data[R]) T(std::forward<F>(args)...);
        R += 1;
    }
    template<typename... F> void emplace_front(F&&... args) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (L <= 0) {
            abort();
        }
        #endif
        L -= 1;
        new(&data[L]) T(std::forward<F>(args)...);
    }
};
template<typename T, int maxsize> class queue_t {
private:
    T* data;
    int L = 0;
    int R = 0;
public:
    queue_t() : data(new T[maxsize]) {}
    ~queue_t() {
        delete[] data;
    }    
    T& front() {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R <= L) {
            abort();
        }
        #endif
        return data[L];
    }
    int size() const {
        return R - L;
    }
    bool empty() const {
        return L == R;
    }
    void clear() {
        L = 0;
        R = 0;
    }
    void push(const T& value) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R >= maxsize) {
            abort();
        }
        #endif
        data[R] = value;
        R += 1;
    }
    template<typename... F> void emplace(F&&... args) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R >= maxsize) {
            abort();
        }
        #endif
        new(&data[R]) T(std::forward<F>(args)...);
        R += 1;
    }
    void pop() {
        #ifdef __SMZ_RUNTIME_CHECK
        if (R <= L) {
            abort();
        }
        #endif
        L += 1;
    }
};
template<typename T, int maxsize> class vector_t {
private:
    T* data;
    int sz = 0;
public:
    vector_t() : data(new T[maxsize]) {}
    ~vector_t() {
        delete[] data;
    }    
    int size() const {
        return sz;
    }
    bool empty() const {
        return sz == 0;
    }
    void clear() {
        sz = 0;
    }
    void push_back(const T& value) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (sz >= maxsize) {
            abort();
        }
        #endif
        data[sz] = value;
        sz += 1;
    }
    template<typename... F> void emplace_back(F&&... args) {
        #ifdef __SMZ_RUNTIME_CHECK
        if (sz >= maxsize) {
            abort();
        }
        #endif
        new(&data[sz]) T(std::forward<F>(args)...);
        sz += 1;
    }
    void pop_back() {
        #ifdef __SMZ_RUNTIME_CHECK
        if (sz <= 0) {
            abort();
        }
        #endif
        sz -= 1;
    }
    T* begin() {
        return data;
    }
    T* end() {
        return data + sz;
    }
};
double runtime() {
    auto now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - start_time);
    return duration.count() / 1e6;
}
namespace io {
    constexpr int MAXBUFFER = 1024 * 1024 * 8;
    char ibuffer[MAXBUFFER], * iptr, obuffer[MAXBUFFER], * optr;
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
        char tmp[32], * now = tmp + 20;
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
    int p[MAXN];
    std::pair<int, int> nodes[MAXM];
    std::vector<std::tuple<int, int, int, int, int, int, std::vector<int>>> business;
    void run() {
        io::start_reading();
        ::n = io::read_int();
        ::m = io::read_int();
        #ifdef __SMZ_RUNTIME_CHECK
        if (n < 2 || n > 200 || m < 1 || m > 1000) {
            abort();
        }
        #endif
        io::start_reading();
        for (int i = 1; i <= n; ++i) {
            p[i] = io::read_int();
            #ifdef __SMZ_RUNTIME_CHECK
            if (p[i] < 0 || p[i] > 20) {
                abort();
            }
            #endif
        }
        for (int i = 1; i <= m; ++i) {
            io::start_reading();
            int x = io::read_int();
            int y = io::read_int();
            nodes[i] = { x, y };
            #ifdef __SMZ_RUNTIME_CHECK
            if (x < 1 || x > n || y < 1 || y > n || x == y) {
                abort();
            }
            #endif
        }
        io::start_reading();
        ::q = io::read_int();
        for (int j = 1; j <= q; ++j) {
            io::start_reading();
            int s = io::read_int();
            int t = io::read_int();
            int length = io::read_int();
            int L = io::read_int();
            int R = io::read_int();
            int V = io::read_int();
            #ifdef __SMZ_RUNTIME_CHECK
            if (s < 1 || s > n || t < 1 || t > n || length < 0 || L < 0 || R > 40 || L > R || V < 0) {
                abort();
            }
            #endif
            io::start_reading();
            std::vector<int> E;
            for (int i = 1; i <= length; ++i) {
                int e = io::read_int();
                E.push_back(e);
                #ifdef __SMZ_RUNTIME_CHECK
                if (e < 1 || e > m) {
                    abort();
                }
                #endif
            }
            business.emplace_back(s, t, length, L, R, V, std::move(E));
        }
    }
    void start() {
        #ifdef __SMZ_RUNTIME_CHECK
        static int64_t id = 0;
        if (id) {
            for (int i = 1; i <= q; ++i) {
                query[i].undo();
            }
            for (int i = 1; i <= n; ++i) {
                if (p[i] != ::p[i]) {
                    abort();
                }
            }
            for (int i = 1; i <= m; ++i) {
                if (edges[i].channel) {
                    abort();
                }
            }
        }
        id += 1;
        #endif
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
            ::query[i].backup.clear();
            for (auto e : E) {
                ::query[i].path.emplace_back(e, L);
            }
            ::query[i].init(::query[i].path);
        }
    }
}
namespace search {
    int first_vis[MAXN];
    int64_t last[MAXQ], visit[MAXN][MAXK], timestamp = 1;
    int dist[MAXN][MAXK], same[MAXN][MAXK];
    std::tuple<int, int, int> father[MAXN][MAXK];
    std::bitset<MAXN> state[MAXN][MAXK];
    deque_t<int, MAXN * MAXK * 4> A, B1, B2, C;
    queue_t<int, MAXN> Q;
    int baseline[MAXN][MAXN];
    inline void preprocess(std::vector<int> deleted) {
        std::sort(deleted.begin(), deleted.end());
        auto iter = std::unique(deleted.begin(), deleted.end());
        deleted.erase(iter, deleted.end());
        for (auto i : deleted) {
            const int start = query[i].to;
            Q.clear();
            for (int i = 1; i <= n; ++i) {
                baseline[start][i] = INF;
            }
            baseline[start][start] = 0;
            Q.push(start);
            while (Q.size()) {
                auto x = Q.front(); Q.pop();
                for (auto [y, _] : G[x]) {
                    if (baseline[start][y] == INF) {
                        baseline[start][y] = baseline[start][x] + 1;
                        Q.push(y);
                    }
                }
            }
        }
    }
    path_t astar(const query_t& qry, int limit=INF) noexcept {
        if (baseline[qry.to][qry.from] == INF) {
            return {};
        }
        timestamp += 2;
        A.clear(); B1.clear(); B2.clear(); C.clear();
        for (int i = 1; i <= n; ++i) {
            first_vis[i] = 0;
        }
        for (int j = 1; j + qry.span <= k; ++j) {
            visit[qry.from][j] = timestamp;
            dist[qry.from][j] = 0;
            state[qry.from][j].reset();
            state[qry.from][j].set(qry.from);
            A.push_back(qry.from | (j << 12));
        }
        int channel = -1;
        int estimate = baseline[qry.to][qry.from];
        while (A.size() || B1.size() || B2.size() || C.size()) {
            while (A.empty()) {
                estimate += 1;
                if (estimate > limit) {
                    return {};
                }
                A.clear();
                //todo：更新合并方式
                while (B1.size() && B2.size()) {
                    if (B1.front() > B2.front()) {
                        A.push_front(B1.front());
                        B1.pop_front();
                    }
                    else {
                        A.push_front(B2.front());
                        B2.pop_front();
                    }
                }
                while (B1.size()) {
                    A.push_front(B1.front());
                    B1.pop_front();
                }
                while (B2.size()) {
                    A.push_front(B2.front());
                    B2.pop_front();
                }
                B1.clear();
                B2.clear();
                std::swap(B2, C);
            }
            auto tmp = A.back(); A.pop_back();
            int x = tmp & 0xFFF, i = tmp >> 12;
            if (visit[x][i] > timestamp) {
                continue;
            }
            if (x == qry.to) [[unlikely]] {
                channel = i;
                break;
            }
            if (p[x] > 0 && !first_vis[x]) {
                C.push_front(x | (0 << 12));
                first_vis[x] = i;
            }
            if (i == 0) {
                const int i = first_vis[x];
                for (int j = 1; j + qry.span <= k; ++j) {
                    if (visit[x][j] < timestamp || dist[x][j] > dist[x][i] + 2) {
                        visit[x][j] = timestamp;
                        dist[x][j] = dist[x][i] + 2;
                        state[x][j] = state[x][i];
                        father[x][j] = {x, i, -1};
                    }
                    A.push_back(x | (j << 12));
                }
                continue;
            }
            int base = dist[x][i] + baseline[qry.to][x];
            visit[x][i] = timestamp + 1;
            const uint64_t mask = (1ull << (i + qry.span + 1)) - (1ull << i);
            for (auto [y, info] : G[x]) {
                if (!info->empty(mask)) {
                    continue;
                }
                if (state[x][i].test(y)) [[unlikely]] {
                    continue;
                }
                if (visit[y][i] < timestamp) {
                    visit[y][i] = timestamp;
                    dist[y][i] = INF;
                }
                const int weight = last[info->index] == timestamp;
                if (dist[y][i] > dist[x][i] + 1) {
                    dist[y][i] = dist[x][i] + 1;
                    state[y][i] = state[x][i];
                    state[y][i].set(y);
                    father[y][i] = {x, std::get<0>(father[x][i]) == x ? std::get<1>(father[x][i]) : i, info->index};
                    int estimate = dist[y][i] + baseline[qry.to][y];
                    if (estimate == base) {
                        A.push_back(y | (i << 12));
                    }
                    else if (estimate == base + 1) {
                        B1.push_back(y | (i << 12));
                    }
                    else {
                        #ifdef __SMZ_RUNTIME_CHECK
                        if (estimate != base + 2) {
                            abort();
                        }
                        #endif
                        C.push_back(y | (i << 12));
                    }
                }
            }
        }
        if (channel == -1) {
            return {};
        }
        const int length = dist[qry.to][channel];
        path_t path;
        path.reserve(length);
        int tmp = channel;
        int node = qry.to;
        while (node != qry.from) {
            auto [prev_node, prev_channel, e] = father[node][channel];
            path.emplace_back(e, channel);
            node = prev_node;
            channel = prev_channel;
            #ifdef __SMZ_RUNTIME_CHECK
            if (node <= 0 || p[node] < -1 || node > n || channel <= 0 || channel > k) {
                abort();
            }
            #endif
        }
        #ifdef __SMZ_RUNTIME_CHECK
        if (path.size() > length) {
            abort();
        }
        #endif
        std::reverse(path.begin(), path.end());
        return path;
    }
    path_t bfs(const query_t& qry, int start_c=1) {
        timestamp += 2;
        int channel = -1;
        static queue_t<std::pair<int, int>, MAXN * MAXK * 2> queue; queue.clear();
        // queue.emplace(qry.from, start_c);
        // visit[qry.from][start_c] = timestamp;
        // while (channel == -1 && !queue.empty()) {
        //     auto [x, i] = queue.front();
        //     queue.pop();
        //     const uint64_t mask = (1ull << (i + qry.span + 1)) - (1ull << i);
        //     for (auto [y, info] : G[x]) {
        //         if (!info->empty(mask)) {
        //             continue;
        //         }
        //         if (visit[y][i] != timestamp) {
        //             visit[y][i] = timestamp;
        //             queue.emplace(y, i);
        //             father[y][i] = {x, i, info->index};
        //             if (y == qry.to) {
        //                 channel = i;
        //                 break;
        //             }
        //         }
        //     }
        // }

        // if (channel == -1) {
        //     std::vector<int> vec;
        //     vec.reserve(k - qry.span);
        //     for(int j = 1; j <= k - qry.span; ++j) {//if (j != start_c) {
        //         vec.push_back(j);
        //     }
        //     // shuffle(vec.begin(), vec.end(), engine);
        //     for (auto j : vec) {
        //         queue.clear();
        //         queue.emplace(qry.from, j);
        //         visit[qry.from][j] = timestamp;
        //         while (channel == -1 && !queue.empty()) {
        //             auto [x, i] = queue.front();
        //             queue.pop();
        //             const uint64_t mask = (1ull << (i + qry.span + 1)) - (1ull << i);
        //             for (auto [y, info] : G[x]) {
        //                 if (!info->empty(mask)) {
        //                     continue;
        //                 }
        //                 if (visit[y][i] != timestamp) {
        //                     visit[y][i] = timestamp;
        //                     queue.emplace(y, i);
        //                     father[y][i] = {x, i, info->index};
        //                     if (y == qry.to) {
        //                         channel = i;
        //                         break;
        //                     }
        //                 }
        //             }
        //         }
        //         if (channel != -1) {
        //             break;
        //         }
        //     }
        // }

        queue.clear();
        for(int j = 1; j <= k - qry.span; ++j) if (j != start_c) {
            queue.emplace(qry.from, j);
            visit[qry.from][j] = timestamp;
        }
        
        while (channel == -1 && !queue.empty()) {
            auto [x, i] = queue.front();
            queue.pop();
            const uint64_t mask = (1ull << (i + qry.span + 1)) - (1ull << i);
            for (auto [y, info] : G[x]) {
                if (!info->empty(mask)) {
                    continue;
                }
                if (visit[y][i] != timestamp) {
                    visit[y][i] = timestamp;
                    queue.emplace(y, i);
                    father[y][i] = {x, i, info->index};
                    if (y == qry.to) {
                        channel = i;
                        break;
                    }
                }
            }
        }

        if (channel == -1) {
            return {};
        }
        path_t path;
        int node = qry.to;
        while (node != qry.from) {
            auto [prev_node, prev_channel, e] = father[node][channel];
            path.emplace_back(e, channel);
            node = prev_node;
            channel = prev_channel;
            #ifdef __SMZ_RUNTIME_CHECK
            if (node <= 0 || p[node] < -1 || node > n || channel <= 0 || channel > k) {
                abort();
            }
            #endif
        }
        std::reverse(path.begin(), path.end());
        return path;
    }
}
template<bool once=true, bool is_baseline=false> transaction_t solve(int e) {
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
    auto iter = std::remove_if(deleted.begin(), deleted.end(), [](int i) {
        return query[i].dead;
    });
    deleted.erase(iter, deleted.end());
    if constexpr (is_baseline) {
        std::sort(deleted.begin(), deleted.end(), [](int x, int y) {
            return query[x].value > query[y].value;
        });
    }
    else {
        std::sort(deleted.begin(), deleted.end(), [](int x, int y) {
            if (query[x].value != query[y].value)
                return query[x].value > query[y].value;
            return query[x].index > query[y].index;
        });
        search::preprocess(deleted);
    }
    std::pair<int64_t, int64_t> best{std::numeric_limits<int64_t>::max(), std::numeric_limits<int64_t>::max()};
    std::vector<std::pair<int, path_t>> answer;
    std::vector<int> order;
    auto proc = [&](const std::vector<int>& indices) {
        int64_t loss = 0, length = 0;
        std::vector<int> updated;
        updated.reserve(indices.size());
        for (auto i : indices) {
            auto c = query[i].path.back().second;
            auto len = query[i].path.size();
            if constexpr (!is_baseline) query[i].undo();
            ::iterations += 1;
            path_t new_path;
            if constexpr (is_baseline) {
                new_path = search::bfs(query[i], c);
            }else{
                new_path = search::astar(query[i], len + 10);
            }
            
            #ifdef __SMZ_RUNTIME_CHECK
            int node = query[i].from;
            std::vector<int> nodes(1, node);
            for (int i = 0; i < new_path.size(); ++i) {
                const auto [e, _] = new_path[i];
                if (p[node] < -1) {
                    abort();
                }
                if (node != edges[e].first && node != edges[e].second) {
                    abort();
                }
                node = (node != edges[e].first ? edges[e].first : edges[e].second);
                for (auto x : nodes) {
                    if (x == node) {
                        abort();
                    }
                }
                nodes.push_back(node);
                for (int j = 0; j < i; ++j) {
                    if (new_path[i].first == new_path[j].first) {
                        abort();
                    }
                }
            }
            if (new_path.size() && node != query[i].to) {
                abort();
            }
            if (query[i].dead) {
                abort();
            }
            #endif
            if (new_path.empty()) {
                loss += query[i].value;
                if constexpr (!is_baseline) query[i].redo();
            }
            else {
                if constexpr (is_baseline) query[i].undo();
                length += new_path.size() * (query[i].span + 1);
                query[i].confirm(std::move(new_path));
                updated.push_back(i);
            }
            if (std::make_pair(loss, length) >= best) {
                break;
            }
        }
        std::vector<std::pair<int, path_t>> result;
        result.reserve(updated.size());
        for (auto i : updated) {
            auto new_path = query[i].cancel();
            result.emplace_back(i, std::move(new_path));
        }
        std::pair<int64_t, int64_t> now{loss, length};
        if (now < best) {
            best = now;
            answer = std::move(result);
            order = indices;
        }
    };
    const double base = runtime();
    const double time_limit = base + (MAXTIME - base) / num_operations;
    if constexpr (once) {
        proc(deleted);
    }
    else {
        if (deleted.size() <= 3) {
            std::vector<int> permutation;
            for (int i = 0; i < deleted.size(); ++i) {
                permutation.push_back(i);
            }
            do {
                std::vector<int> indices;
                for (auto i : permutation) {
                    indices.push_back(deleted[i]);
                }
                proc(indices);
            } while (runtime() < time_limit && std::next_permutation(permutation.begin(), permutation.end()));
        }
        else {
            proc(deleted);
            std::bernoulli_distribution bernoulli(std::pow(deleted.size(), -1.0 / 3));
            while (runtime() < time_limit) {
                std::vector<int> indices = order;
                for (int i = 1; i < indices.size(); ++i) if (bernoulli(engine)) {
                    int j = std::rand() % i;
                    std::swap(indices[i], indices[j]);
                }
                proc(indices);
            }
        }
    }
    static uint64_t flag[MAXQ], timestamp = 1;
    timestamp += 1;
    std::vector<int> dead;
    for (const auto &[i, _] : answer) {
        flag[i] = timestamp;
    }
    for (auto i : deleted) {
        if (flag[i] != timestamp) {
            dead.push_back(i);
        }
    }
    return transaction_t{best.first, std::move(answer), std::move(dead)};
}

namespace union_set {
    int fa[MAXN];
    void init(){
        for(int i = 1; i <= n; i++){
            fa[i] = i;
        }
    }

    int find(int x){
        if(fa[x] != x) fa[x] = find(fa[x]);
        return fa[x];
    }

    std::vector <int> gen(int64_t *visit, int64_t timestamp){
        std::vector <int> indices(m);
        std::vector <int> ret;
        for(int i = 1; i <= m; i++){
            indices[i - 1] = i;
        }
        shuffle(indices.begin(), indices.end(), engine);

        init();
        for(auto i : indices){
            if (visit[i] == timestamp) continue;
            int x = edges[i].first, y = edges[i].second;
            if(find(x) == find(y)){
                ret.push_back(i);
            }else{
                fa[find(x)] = find(y);
            }
        }
        return ret;
    }
}

void generate() { //输出瓶颈断边场景的交互部分
    static int64_t visit[MAXM];
    static int64_t timestamp = 1;
    auto check = [](const auto& deleted) {
        auto jaccard = [](const auto& A, const auto& B) {
            double x = 0, y = A.size();
            timestamp += 1;
            for (auto [_, v] : A) {
                visit[v] = timestamp;
            }
            for (auto [_, v] : B) {
                if (visit[v] == timestamp) {
                    x += 1;
                }
                else {
                    y += 1;
                }
            }
            return x / y;
        };
        for (int i = 0; i < pretests.size(); ++i) {
            if (jaccard(pretests[i], deleted) > MAXJACCARD) {
                return false;
            }
        }
        return true;
    };
    testcase::start();
    std::vector<int> indices;
    for (int i = 1; i <= m; ++i) {
        indices.push_back(i);
    }
    std::uniform_real_distribution<double> eps(0, 1e-5);
    std::priority_queue<item_t> queue;
    case_t init;
    for (int i = 1; i <= m; ++i) {
        init.emplace_back(0.0, i);
    }
    std::shuffle(init.begin(), init.end(), engine);
    if (init.size() > MAXC) {
        init.resize(MAXC);
    }
    queue.emplace(0.0, DEGREE, std::move(init));
    double last = runtime();
    std::vector<std::pair<double, case_t>> cases;
    vector_t<int, MAXC> position[MAXC];
    for (;;) {
        double now = runtime();
        if (now > MAXGENTIME) {
            break;
        }
        std::vector<std::pair<double, int>> deleted;
        timestamp += 1;
        while (queue.size() > 1 && queue.top().counter <= 0) {
            auto item = queue.top();
            queue.pop();
            cases.emplace_back(item.value, std::move(item.sequence));
        }
        const auto& [_, cnt, best] = queue.top();
        const_cast<int&>(cnt) -= 1;
        deleted = best;
        const int maxsize = std::min<int>(deleted.size() + MAXC * EXPAND_RATIO, MAXC);
        for (int i = 0; i < deleted.size(); ++i) {
            deleted[i].first += eps(engine);
        }
        std::vector<int> order(deleted.size());
        for (int i = 0; i < deleted.size(); ++i) {
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [&](int x, int y) {
            return deleted[x].first < deleted[y].first;
        });
        int r = std::max(std::pow(deleted.size(), TRANSFORM_EXP), 1.0);
        for (int i = 0; i < r; ++i) {
            deleted[order[i]].second = -1;
        }
        auto iter = std::remove_if(deleted.begin(), deleted.end(), [](auto pair) {
            return pair.second == -1;
        });
        deleted.erase(iter, deleted.end());
        // std::shuffle(indices.begin(), indices.end(), engine);
        for (auto [v, i] : deleted) {
            visit[i] = timestamp;
        }
        auto indices = union_set::gen(visit, timestamp);
        for (auto i : indices) if (visit[i] != timestamp) {
            deleted.emplace_back(0, i);
            if (deleted.size() >= MAXC) {
                break;
            }
        }
        //constexpr double probability = (1 - MAXJACCARD) / (1 + MAXJACCARD);
        std::bernoulli_distribution bernoulli(1.0 / r);
        std::uniform_int_distribution<int> uniform(0, MAXC - 1);
        for (int i = 0; i < MAXC; ++i) {
            position[i].clear();
        }
        for (int i = (int)deleted.size() - 1; i >= 0; --i) {
            int j = i;
            if (bernoulli(engine)) {
                j = uniform(engine);
            }
            position[j].push_back(deleted[i].second);
        }
        deleted.clear();
        for (int i = 0; i < MAXC; ++i) {
            for (auto j : position[i]) {
                deleted.emplace_back(0.0, j);
            }
        }
        if (deleted.size() > maxsize) {
            deleted.resize(maxsize);
        }
        testcase::start();
        for (int j = 0; j < deleted.size(); ++j) {
            int e = deleted[j].second;
            auto transaction = solve<true, false>(e);
            deleted[j].first = -transaction.loss;
            transaction.commit();
        }
        testcase::start();
        for (int j = 0; j < deleted.size(); ++j) {
            int e = deleted[j].second;
            auto transaction = solve<true, true>(e);
            deleted[j].first += transaction.loss;
            transaction.commit();
        }
        int idx = 0;
        double delta = 0, mx = 0;
        for (int j = 0; j < deleted.size(); ++j) {
            delta += deleted[j].first;
            if (delta > mx) { //todo：>=
                mx = delta;
                idx = j;
            }
        }
        deleted.resize(idx + 1);
        queue.emplace(mx, DEGREE, std::move(deleted));
    }
    while (queue.size()) {
        auto item = queue.top();
        queue.pop();
        cases.emplace_back(item.value, std::move(item.sequence));
    }
    #ifdef __SMZ_RUNTIME_CHECK
    for (auto& [_, deleted] : cases) {
        timestamp += 1;
        for (auto [_, i] : deleted) {
            if (i <= 0 || i > m) {
                abort();
            }
            if (visit[i] == timestamp) {
                abort();
            }
            visit[i] = timestamp;
        }
    }
    #endif
    std::sort(cases.begin(), cases.end(), [](const auto& x, const auto& y) {
        return x.first > y.first;
    });
    int used = -1;
    for (int i = 0; i < cases.size(); ++i) {
        auto& [_, deleted] = cases[i];
        if (check(deleted)) {
            pretests.push_back(std::move(deleted));
            if (pretests.size() == MAXT1) {
                used = i;
                break;
            }
        }
    }
    io::start_writing();
    io::write_int(pretests.size());
    io::flush();
    for(const auto& deleted : pretests){
        io::start_writing();
        io::write_int(deleted.size());
        io::newline();
        for (auto [_, e] : deleted) {
            io::write_int(e);
        }
        io::flush();
    }
    #ifdef __SMZ_NATIVE_TEST
    if (pretests.size() > MAXT1) {
        abort();
    }
    for (const auto& deleted : pretests) {
        if (deleted.size() > MAXC) {
            abort();
        }
    }
    print("Data Generated.");
    print("Cases size: ", used + 1, "/", cases.size());
    print("Max: ", cases[0].first);
    double sum = 0, total = 0;
    for (int j = 1; j <= q; ++j) {
        total += query[j].value;
    }
    for (const auto& deleted : pretests) {
        for (auto [v, _] : deleted) {
            sum += v;
        }
    }
    print("Score different: ", sum * 10000.0 / total);
    #endif
}
int main() { // 244664 388723 44496.7(1496772)
#ifdef __SMZ_NATIVE_TEST
    std::ignore = freopen("cases/testcase8.in", "r", stdin);
    std::ignore = freopen("output.txt", "w", stdout);
#endif
    testcase::run();
    generate();
    //输出瓶颈断边场景的交互部分
    io::start_reading();
    int T = io::read_int();
    int idx = 0;
    #ifdef __SMZ_NATIVE_TEST
    T += pretests.size();
    if (pretests.empty()) {
        abort();
    }
    #endif
    double score = 0;
    while (T) {
        testcase::start();
        std::vector<std::pair<double, int>> data;
#ifdef __SMZ_NATIVE_TEST
        uint64_t total = 0, rest = 0;
        for (int i = 1; i <= q; ++i) {
            total += query[i].value;
        }
        if (idx < pretests.size()) {
            data = pretests[idx];
            data.push_back({0.0, -1});
            std::reverse(data.begin(), data.end());
        }
#endif
        const int maxfail = std::min(m, MAXC);
        if (num_operations > T * maxfail) {
            num_operations = T * maxfail;
        }
        for (;;) {
            int e;
            if (data.size()) {
                e = data.back().second;
                data.pop_back();
            }
            else {
                io::start_reading();
                e = io::read_int();
            }
            if (e == -1) {
                break;
            }
            auto transaction = solve(e);
            io::start_writing();
            io::write_int((int)transaction.result.size());
            io::newline();
            for (const auto& [i, path] : transaction.result) {
                io::write_int(i);
                io::write_int((int)path.size());
                io::newline();
                for (auto [e, c] : path) {
                    io::write_int(e);
                    io::write_int(c);
                    io::write_int(c + query[i].span);
                }
                io::newline();
            }
            io::flush();
            transaction.commit();
            edges[e].deleted = true;
            num_operations -= 1;
        }
        T -= 1;
#ifdef __SMZ_NATIVE_TEST
        for (int i = 1; i <= q; ++i) if (!query[i].dead) {
            rest += query[i].value;
        }
        score += rest * 10000.0 / total;
        idx += 1;
#endif
    }
#ifdef __SMZ_NATIVE_TEST
    print("Score: ", (int)score);
    print("Runtime: ", runtime());
    print("Iterations: ", iterations);
#endif
    return 0;
}
