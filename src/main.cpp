#include <bits/stdc++.h>
constexpr int K = 44;
constexpr int MAXN = 256;
constexpr int MAXM = 1024;
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
    inline void flush() {
        if (optr != obuffer) {
            optr[-1] = '\n';
        }
        fwrite(obuffer, 1, optr - obuffer, stdout);
        fflush(stdout);
    }
}
int N, M, J, P[MAXN];
int main() {
    #ifdef __SMZ_RUNTIME_CHECK
    std::ignore = freopen("../release/testcase1.in", "r", stdin);
    #endif
    io::start_reading();
    ::N = io::read_int();
    ::M = io::read_int();
    io::start_reading();
    for (int i = 1; i <= N; ++i) {
        P[i] = io::read_int();
    }
    for (int i = 1; i <= M; ++i) {
        io::start_reading();
        int x = io::read_int();
        int y = io::read_int();
    }
    io::start_reading();
    ::J = io::read_int();
    for (int j = 1; j <= J; ++j) {
        io::start_reading();
        int s = io::read_int();
        int t = io::read_int();
        int m = io::read_int();
        int L = io::read_int();
        int R = io::read_int();
        int V = io::read_int();
        io::start_reading();
        for (int i = 1; i <= m; ++i) {
            int e = io::read_int();
        }
    }
    io::start_reading();
    int T = io::read_int();
    while (T--) {
        for (;;) {
            io::start_reading();
            int e = io::read_int();
            if (e == -1) {
                break;
            }
            io::start_writing();
            io::write_int(0);
            io::flush();
        }
    }
    return 0;
}
