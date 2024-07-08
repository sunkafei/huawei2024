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
constexpr int INF = 1 << 30;
constexpr int k = 40;
constexpr int MAXK = 44;
constexpr int MAXN = 256;
constexpr int MAXM = 1024;
constexpr int MAXQ = 6000;
constexpr int MAXTIME = 89;
int64_t iterations = 0;
int num_operations = INF;
int n, m, q, p[MAXN];
using path_t = std::vector<std::pair<int, int>>;
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
std::mt19937 engine;
std::vector<std::pair<int, edge_t*>> G[MAXN];
const auto start_time = std::chrono::steady_clock::now();
template<typename... T> void print(const T&... sth) {
#ifdef __SMZ_NATIVE_TEST
    (..., (std::cerr << sth << " ")) << std::endl;
#endif
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
        L -= 1;
        data[L] = value;
    }
    void push_back(const T& value) {
        data[R] = value;
        R += 1;
    }
    void pop_front() {
        L += 1;
    }
    void pop_back() {
        R -= 1;
    }
    template<typename... F> void emplace_back(F&&... args) {
        new(&data[R]) T(std::forward<F>(args)...);
        R += 1;
    }
    template<typename... F> void emplace_front(F&&... args) {
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
        data[R] = value;
        R += 1;
    }
    template<typename... F> void emplace(F&&... args) {
        new(&data[R]) T(std::forward<F>(args)...);
        R += 1;
    }
    void pop() {
        L += 1;
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
    int64_t last[MAXQ], visit[MAXN][MAXK], timestamp = 1;
    int dist[MAXN][MAXK], same[MAXN][MAXK];
    int first_vis[MAXN];
    std::tuple<int, int, int> father[MAXN][MAXK];
    deque_t<std::pair<int, int>, MAXN * MAXK> queue;
    template<typename T> void preprocess(const T&) {
        
    }
    inline path_t search(const query_t& qry, const path_t& prev) noexcept {
        static std::bitset<MAXN> state[MAXN][MAXK];
        timestamp += 2;
        queue.clear();
        for (auto [e, _] : prev) {
            last[e] = timestamp;
        }
        for (int i = 1; i <= n; ++i) {
            first_vis[i] = false;
        }
        for (int j = k - qry.span; j > 0; --j) {
            visit[qry.from][j] = timestamp;
            same[qry.from][j] = 0;
            dist[qry.from][j] = 0;
            state[qry.from][j].reset();
            state[qry.from][j].set(qry.from);
            queue.emplace_back(qry.from, j);
        }
        int channel = -1;
        while (!queue.empty()) {
            auto [x, i] = queue.front();
            queue.pop_front();
            if (x == qry.to) [[unlikely]] {
                channel = i;
                break;
            }
            if (visit[x][i] > timestamp) continue;
            if (p[x] > 0 && !first_vis[x]) {
                first_vis[x] = true;
                for (int j = 1; j + qry.span <= k; ++j) {
                    if (visit[x][j] < timestamp || dist[x][j] > dist[x][i]) {
                        visit[x][j] = timestamp;
                        same[x][j] = same[x][i];
                        dist[x][j] = dist[x][i];
                        state[x][j] = state[x][i];
                        father[x][j] = {x, i, -1};
                    }
                    queue.emplace_front(x, j);
                }
                continue;
            }
            visit[x][i] = timestamp + 1;
            const uint64_t mask = (1ull << (i + qry.span + 1)) - (1ull << i);
            for (auto [y, info] : G[x]) {
                if (!info->empty(mask)) {
                    continue;
                }
                if (state[x][i].test(y)) [[unlikely]] {
                    continue;
                }
                const int weight = last[info->index] == timestamp;
                if (visit[y][i] < timestamp || dist[y][i] > dist[x][i] + 1) {
                    visit[y][i] = timestamp;
                    same[y][i] = same[x][i] + weight;
                    dist[y][i] = dist[x][i] + 1;
                    state[y][i] = state[x][i];
                    state[y][i].set(y);
                    queue.emplace_back(y, i);
                    father[y][i] = {x, std::get<0>(father[x][i]) == x ? std::get<1>(father[x][i]) : i, info->index};
                }
                else if (dist[y][i] == dist[x][i] + 1 && same[y][i] < same[x][i] + weight) {
                    same[y][i] = same[x][i] + weight;
                    state[y][i] = state[x][i];
                    state[y][i].set(y);
                    father[y][i] = {x, std::get<0>(father[x][i]) == x ? std::get<1>(father[x][i]) : i, info->index};
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
    std::sort(deleted.begin(), deleted.end(), [](int x, int y) {
        if (query[x].value != query[y].value)
            return query[x].value > query[y].value;
        return query[x].index > query[y].index;
    });
    search::preprocess(deleted);
    std::tuple<int64_t, int64_t> best{std::numeric_limits<int64_t>::max(), std::numeric_limits<int64_t>::max()};
    std::vector<std::pair<int, path_t>> answer;
    std::vector<int> order;
    auto proc = [&](const std::vector<int>& indices) {
        int64_t loss = 0, length = 0;
        std::vector<int> updated;
        updated.reserve(indices.size());
        for (auto i : indices) {
            query[i].undo();
            const auto& prev = query[i].backup;
            ::iterations += 1;
            auto new_path = search::search(query[i], prev);
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
                query[i].redo();
            }
            else {
                length += new_path.size() * (query[i].span + 1);
                query[i].confirm(std::move(new_path));
                updated.push_back(i);
            }
            if (std::make_tuple(loss, length) >= best) {
                break;
            }
        }
        std::vector<std::pair<int, path_t>> result;
        result.reserve(updated.size());
        for (auto i : updated) {
            auto new_path = query[i].cancel();
            result.emplace_back(i, std::move(new_path));
        }
        std::tuple<int64_t, int64_t> now{loss, length};
        if (now < best) {
            best = now;
            answer = std::move(result);
            order = indices;
        }
    };
    const double base = runtime();
    const double time_limit = base + (MAXTIME - base) / num_operations;
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
        std::bernoulli_distribution bernoulli(1.0 / std::sqrt(deleted.size()));
        while (runtime() < time_limit) {
            std::vector<int> indices = order;
            for (int i = 1; i < indices.size(); ++i) if (bernoulli(engine)) {
                int j = std::rand() % i;
                std::swap(indices[i], indices[j]);
            }
            proc(indices);
        }
    }
    static uint64_t flag[MAXQ], timestamp = 1;
    timestamp += 1;
    std::vector<int> ret;
    for (auto iter = answer.begin(); iter != answer.end(); ++iter) {
        auto [i, new_path] = std::move(*iter);
        flag[i] = timestamp;
        ret.push_back(i);
        query[i].replace(std::move(new_path));
    }
    for (auto i : deleted) {
        if (flag[i] != timestamp) {
            query[i].dead = true;
        }
    }
#ifdef __SMZ_RUNTIME_CHECK
    for (auto i : ret) {
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
    return ret;
}
int main() {
#ifdef __SMZ_NATIVE_TEST
    std::ignore = freopen("../release/testcase2.in", "r", stdin);
    std::ignore = freopen("../release/output.txt", "w", stdout);
#endif
    testcase::run();

    //输出瓶颈断边场景的交互部分
    io::start_writing();
    const int T1 = m > 120 ? 50 : 0;
    io::write_int(T1);
    io::flush();
    std::uniform_int_distribution<int> gen(1, m);
    std::mt19937 mt(20140920);
    for (int i = 0; i < T1; ++i) {
        int c = 50;
        io::start_writing();
        io::write_int(c);
        io::newline();
        std::unordered_set<int> visit;
        for (int j = 0; j < c; ++j) {
            int e = gen(mt);
            while (visit.count(e)) {
                e = gen(mt);
            }
            visit.insert(e);
            io::write_int(e);
        }
        io::flush();
    }

    //输出瓶颈断边场景的交互部分
    io::start_reading();
    int T = io::read_int();
    double score = 0;
    while (T) {
        testcase::start();
#ifdef __SMZ_NATIVE_TEST
        uint64_t total = 0, rest = 0;
        for (int i = 1; i <= q; ++i) {
            total += query[i].value;
        }
#endif
        const int maxfail = std::min(m, 50);
        if (num_operations > T * maxfail) {
            num_operations = T * maxfail;
        }
        for (;;) {
            io::start_reading();
            int e = io::read_int();
            if (e == -1) {
                break;
            }
            auto indices = solve(e);
            io::start_writing();
            io::write_int((int)indices.size());
            io::newline();
            for (auto i : indices) {
                io::write_int(i);
                io::write_int((int)query[i].path.size());
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
            num_operations -= 1;
        }
        T -= 1;
#ifdef __SMZ_NATIVE_TEST
        for (int i = 1; i <= q; ++i) if (!query[i].dead) {
            rest += query[i].value;
        }
        score += rest * 10000.0 / total;
#endif
    }
#ifdef __SMZ_NATIVE_TEST
    print("Score: ", (int)score);       //8204340   480152
    print("Runtime: ", runtime());
    print("Iterations: ", iterations);  //53681198  11316254
#endif
    return 0;
}