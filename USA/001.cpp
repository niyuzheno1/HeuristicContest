const double TL2 = 4.9;
const double TL1 = 4.7;
const double TL0 = 2.0;
#include <bits/stdc++.h>
using namespace std;
#define F0(i,n) for (int i=0; i<n; i++)
#define FR(i,n) for (int i=n-1; i>=0; i--)
#define F1(i,n) for (int i=1; i<=n; i++)
#define CL(a,x) memset(x, a, sizeof(x));
#define SZ(x) ((int)x.size())
const int inf = 1000009;
const double pi = acos(-1.0);
typedef pair<int, int> pii;
typedef long long ll;
typedef unsigned long long ull;
const double EPS = 1e-9;
#define PR(x) cerr << #x << "=" << x << endl

inline ll GetTSC() {
    ll lo, hi;
    asm volatile ("rdtsc": "=a"(lo), "=d"(hi));
    return lo + (hi << 32);
}
inline double GetSeconds() {
    return GetTSC() / 2.99e9;
}

const int MAX_RAND = 1 << 30;
struct Rand {
    ll x, y, z, w, o;
    Rand(ll seed) { reseed(seed); o = 0; }
    inline void reseed(ll seed) { x = 0x498b3bc5 ^ seed; y = 0; z = 0; w = 0;  F0(i, 20) mix(); }
    inline void mix() { ll t = x ^ (x << 11); x = y; y = z; z = w; w = w ^ (w >> 19) ^ t ^ (t >> 8); }
    inline ll rand() { mix(); return x & (MAX_RAND - 1); }
    inline int nextInt(int n) { return rand() % n; }
    inline int nextInt(int L, int R) { return rand()%(R - L + 1) + L; }
    inline int nextBool() { if (o < 4) o = rand(); o >>= 2; return o & 1; }
    double nextDouble() { return rand() * 1.0 / MAX_RAND; }
};
Rand my(2019);
double saveTime;
double t_o[20];
ll c_o[20];
void Init() {
    saveTime = GetSeconds();
    F0(i, 20) t_o[i] = 0.0;
    F0(i, 20) c_o[i] = 0;
}
double Elapsed() { return GetSeconds() - saveTime; }
void Report() {
    cerr << "-------------------------------------" << endl;
    cerr << "Elapsed time: " << Elapsed() << " sec" << endl;
    double total = 0.0; F0(i, 20) { if (t_o[i] > 0) cerr << "t_o[" << i << "] = " << t_o[i] << endl; total += t_o[i]; } cerr << endl; //if (total > 0) cerr << "Total spent: " << total << endl;
    F0(i, 20) if (c_o[i] > 0) cerr << "c_o[" << i << "] = " << c_o[i] << endl;
    cerr << "-------------------------------------" << endl;
}
struct AutoTimer {
    int x;
    double t;
    AutoTimer(int x) : x(x) {
        t = Elapsed();
    }
    ~AutoTimer() {
        t_o[x] += Elapsed() - t;
    }
};
#define AT(i) AutoTimer a##i(i)
//#define AT(i)

// CONSTANTS
const int N = 256;
const int M = 10000;
typedef pair<pii, pii> Rect;
const int DX[] = {-1,0,1,0};
const int DY[] = {0,-1,0,1};
int opp(int k) { return k^2; }
Rect rect[N], nrect[N];
double maxd;
int n;
int X[N], Y[N], r[N], lck[N];
double bscore;
Rect bans[N];
double delta, act_score, score;
int bad, dir;
int X1[N], Y1[N], X2[N], Y2[N], s[N];
int SUBX1, SUBX2, SUBY1, SUBY2;
const int LOGN = 1 << 16;
double logs[LOGN];
int Area(Rect r) { return (r.second.first - r.first.first + 1) * (r.second.second - r.first.second + 1); }
int Lenx(Rect r) { return (r.second.first - r.first.first + 1); }
int Leny(Rect r) { return (r.second.second - r.first.second + 1); }
bool Inside(Rect r) { return r.first.first >= 0 && r.first.second >= 0 && r.second.first < M && r.second.second < M; }
Rect Grow(Rect r, int k, int d) {
    if (k == 0) r.first.first -= d;
    if (k == 1) r.first.second -= d;
    if (k == 2) r.second.first += d;
    if (k == 3) r.second.second += d;
    return r;
}
int RoomX(Rect r, int t) {
    int dx = Lenx(r);
    int dy = Leny(r);
    return t / dy - dx;
}
int RoomY(Rect r, int t) {
    int dx = Lenx(r);
    int dy = Leny(r);
    return t / dx - dy;
}
int RoomK(Rect r, int k, int t) {
    if (k&1) return RoomY(r, t);
    return RoomX(r, t);
}
bool MaxedOut(Rect r, int t) {
    int dx = Lenx(r);
    int dy = Leny(r);
    return dx * (dy + 1) > t && dy * (dx + 1) > t;
}
bool Growable(Rect r, int k, int t) {
    r = Grow(r, k, 1);
    return Inside(r) && Area(r) <= t;
}
bool Growable(Rect r, int k) {
    r = Grow(r, k, 1);
    return Inside(r);
}
bool Contains(Rect r, int X, int Y) {
    return r.first.first <= X && X <= r.second.first &&
            r.first.second <= Y && Y <= r.second.second;
}
bool Contains(Rect r1, Rect r2) {
    return r1.first.first <= r2.first.first &&
            r1.first.second <= r2.first.second &&
            r2.second.first <= r1.second.first &&
            r2.second.second <= r1.second.second;
}
bool Contains(Rect r1, int X1, int Y1, int X2, int Y2) {
    return r1.first.first <= X1 &&
            r1.first.second <= Y1 &&
            X2 <= r1.second.first &&
            Y2 <= r1.second.second;
}
bool IntersectsX(Rect r1, Rect r2) {
    if (r1.first.first > r2.second.first || r2.first.first > r1.second.first) return false;
    return true;
}
bool IntersectsY(Rect r1, Rect r2) {
    if (r1.first.second > r2.second.second || r2.first.second > r1.second.second) return false;
    return true;
}
bool Intersects(Rect r1, Rect r2) {
    if (r1.first.first > r2.second.first || r2.first.first > r1.second.first) return false;
    if (r1.first.second > r2.second.second || r2.first.second > r1.second.second) return false;
    return true;
}
void Print(Rect r) { cerr << "(" << r.first.first << "," << r.first.second << ")-(" << r.second.first << "," << r.second.second << ")"; }

double raw(int x, int y) {
    double ratio = (1.0 - 1.0 * min(x, y) / max(x, y));
    return 1 - ratio * ratio;
}

double upper(int x, int y) {
    double ratio = (1.0 - 1.0 * min(x, y) / y);
    return 1 - ratio * ratio;
}

int biases[N];
int biasind[N];
int biased(int k) {
    int bs = 0;
    F0(i, k) {
        bs += biases[i];
    }
    if (!bs) return -1;
    int r = my.nextInt(bs);
    F0(i, k) if (r < biases[i]) return i; else r -= biases[i];
    return -1;
}

int id[N];
int sq(int x) { return x * x; }

double NaiveScore() {
    double ret = 0.0;
    F0(i, n) {
        double rw = 0;
        if (X1[i] <= X[i] && X[i] <= X2[i] && Y1[i] <= Y[i] && Y[i] <= Y2[i]) {
            rw = raw((X2[i]-X1[i]+1) * (Y2[i]-Y1[i]+1), r[i]);
        }
        ret += rw;
    }
    return ret;
}

const int CR = 1024;
const int LG = 10;
const int MAXG = M/CR+10;
int tn[MAXG][MAXG], t[MAXG][MAXG][128];
int BX1, BX2, BY1, BY2;
int bx1, bx2, by1, by2;
int ox1, oy1, ox2, oy2, oa;
int qox1, qoy1, qox2, qoy2;
int fx[M], fy[M], nx[N], ny[N];
int undecided[N], unn;

void PrepG(int i) {
    BX1 = X1[i] >> LG;
    BY1 = Y1[i] >> LG;
    BX2 = X2[i] >> LG;
    BY2 = Y2[i] >> LG;
}

void PrepG() {
    BX1 = bx1 >> LG;
    BY1 = by1 >> LG;
    BX2 = bx2 >> LG;
    BY2 = by2 >> LG;
}

void AddTValue(int x, int y, int i) {
    t[x][y][tn[x][y]++] = i;
}

void RemTValue(int x, int y, int i) {
    F0(k, tn[x][y]) if (t[x][y][k] == i) {
        for (int l = k; l + 1 < tn[x][y]; l++) t[x][y][l] = t[x][y][l + 1];
        tn[x][y]--;
        return;
    }
}

struct small_set {
    int wh[N];
    int n;
    int a[N];
    small_set() {
    }
    void erase(int x) {
        int k = wh[x];
        swap(a[--n], a[k]);
        wh[a[k]] = k;
    }
    void insert(int x) {
        a[n] = x;
        wh[x] = n++;
    }
};

small_set coverx[N], covery[N];
void AddT(int i) {
    PrepG(i);
    for (int X = BX1; X <= BX2; X++) for (int Y = BY1; Y <= BY2; Y++) AddTValue(X, Y, i);
}

void RemT(int i) {
    PrepG(i);
    for (int X = BX1; X <= BX2; X++) for (int Y = BY1; Y <= BY2; Y++) RemTValue(X, Y, i);
}

void ChangeCovers(int i) {
    if (X1[i] != qox1) {
        if (X1[i] < qox1) {
            for (int j = fx[X1[i]]; X[j] < qox1; j = nx[j]) coverx[j].insert(i);
        } else {
            for (int j = fx[qox1]; X[j] < X1[i]; j = nx[j]) coverx[j].erase(i);
        }
    }
    if (X2[i] != qox2) {
        if (X2[i] > qox2) {
            for (int j = fx[qox2+1]; X[j] <= X2[i]; j = nx[j]) coverx[j].insert(i);
        } else {
            for (int j = fx[X2[i]+1]; X[j] <= qox2; j = nx[j]) coverx[j].erase(i);
        }
    }
    if (Y1[i] != qoy1) {
        if (Y1[i] < qoy1) {
            for (int j = fy[Y1[i]]; Y[j] < qoy1; j = ny[j]) covery[j].insert(i);
        } else {
            for (int j = fy[qoy1]; Y[j] < Y1[i]; j = ny[j]) covery[j].erase(i);
        }
    }
    if (Y2[i] != qoy2) {
        if (Y2[i] > qoy2) {
            for (int j = fy[qoy2+1]; Y[j] <= Y2[i]; j = ny[j]) covery[j].insert(i);
        } else {
            for (int j = fy[Y2[i]+1]; Y[j] <= qoy2; j = ny[j]) covery[j].erase(i);
        }
    }
}

int seenZ[N], Z;
void GenBoundaryFast(int i) {
    bx1 = SUBX1; by1 = SUBY1;
    bx2 = SUBX2; by2 = SUBY2;

    unn = 0;
    Z++;
    F0(u, coverx[i].n) {
        int j = coverx[i].a[u];
        seenZ[j] = Z;
        if (Y2[j] < Y[i]) by1 = max(by1, Y2[j]+1);
        else by2 = min(by2, Y1[j]-1);
    }
    F0(u, covery[i].n) {
        int j = covery[i].a[u];
        seenZ[j] = Z;
        if (X2[j] < X[i]) bx1 = max(bx1, X2[j]+1);
        else bx2 = min(bx2, X1[j]-1);
    }

    PrepG();
    for (int x = BX1; x <= BX2; x++) for (int y = BY1; y <= BY2; y++) {
        F0(k, tn[x][y]) {
            int j = t[x][y][k];
            if (seenZ[j] != Z) {
                seenZ[j] = Z;
                if (X1[j] <= X2[i] && X2[j] >= X1[i]) {
                    if (Y2[j] < Y[i]) by1 = max(by1, Y2[j]+1);
                    else by2 = min(by2, Y1[j]-1);
                }
                else if (Y1[j] <= Y2[i] && Y2[j] >= Y1[i]) {
                    if (X2[j] < X[i]) bx1 = max(bx1, X2[j]+1);
                    else bx2 = min(bx2, X1[j]-1);
                } else {
                    undecided[unn++] = j;
                }
            }
        }
    }

    int dx = 0, dy = 0;
    F0(uuu, unn) {
        int j = undecided[uuu];
        if (X1[j] > X2[i] && (dx = bx2-X1[j]+1) > 0) {
            if (Y1[j] > Y2[i] && (dy = by2-Y1[j]+1) > 0) {
                if (my.nextInt(dx+dy) < dx) {
                    by2 -= dy;
                } else {
                    bx2 -= dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                if (my.nextInt(dx+dy) < dx) {
                    by1 += dy;
                } else {
                    bx2 -= dx;
                }
            }
        }
        if (X1[j] < X2[i] && (dx = X2[j]+1-bx1)>0) {
            if (Y1[j] > Y2[i] && (dy = by2-Y1[j]+1) > 0) {
                if (my.nextInt(dx+dy) < dx) {
                    by2 -= dy;
                } else {
                    bx1 += dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                if (my.nextInt(dx+dy) < dx) {
                    by1 += dy;
                } else {
                    bx1 += dx;
                }
            }
        }
    }
}

void GenBoundaryFaster(int i) {
    bx1 = SUBX1; by1 = SUBY1;
    bx2 = SUBX2; by2 = SUBY2;

    unn = 0;
    Z++;
    F0(u, coverx[i].n) {
        int j = coverx[i].a[u];
        seenZ[j] = Z;
        if (Y2[j] < Y[i]) by1 = max(by1, Y2[j]+1);
        else by2 = min(by2, Y1[j]-1);
    }
    F0(u, covery[i].n) {
        int j = covery[i].a[u];
        seenZ[j] = Z;
        if (X2[j] < X[i]) bx1 = max(bx1, X2[j]+1);
        else bx2 = min(bx2, X1[j]-1);
    }

    PrepG();
    for (int x = BX1; x <= BX2; x++) for (int y = BY1; y <= BY2; y++) {
        F0(k, tn[x][y]) {
            int j = t[x][y][k];
            if (seenZ[j] != Z) {
                seenZ[j] = Z;
                undecided[unn++] = j;
            }
        }
    }

    int dx = 0, dy = 0;
    F0(uuu, unn) {
        int j = undecided[uuu];
        if (X1[j] > X2[i] && (dx = bx2-X1[j]+1) > 0) {
            if (Y1[j] > Y2[i] && (dy = by2-Y1[j]+1) > 0) {
                biases[0] = dx * (by2 - by1 + 1);
                biases[1] = dy * (bx2 - bx1 + 1);
                if (my.nextInt(biases[0]+biases[1]) < biases[0]) {
                    by2 -= dy;
                } else {
                    bx2 -= dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                biases[0] = dx * (by2 - by1 + 1);
                biases[1] = dy * (bx2 - bx1 + 1);
                if (my.nextInt(biases[0]+biases[1]) < biases[0]) {
                    by1 += dy;
                } else {
                    bx2 -= dx;
                }
            }
        }
        if (X1[j] < X2[i] && (dx = X2[j]+1-bx1)>0) {
            if (Y1[j] > Y2[i] && (dy = by2-Y1[j]+1) > 0) {
                biases[0] = dx * (by2 - by1 + 1);
                biases[1] = dy * (bx2 - bx1 + 1);
                if (my.nextInt(biases[0]+biases[1]) < biases[0]) {
                    by2 -= dy;
                } else {
                    bx1 += dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                biases[0] = dx * (by2 - by1 + 1);
                biases[1] = dy * (bx2 - bx1 + 1);
                if (my.nextInt(biases[0]+biases[1]) < biases[0]) {
                    by1 += dy;
                } else {
                    bx1 += dx;
                }
            }
        }
    }
}

void GenBoundary2(int i, int cnt) {
    Rect bb;
    F0(k, cnt) {
        GenBoundaryFaster(i);
        Rect cc(pii(bx1, by1), pii(bx2, by2));
        if (k == 0 || Area(bb) < Area(cc)) bb = cc;
    }
    bx1 = bb.first.first;
    by1 = bb.first.second;
    bx2 = bb.second.first;
    by2 = bb.second.second;
}

void UpdateBest() {
    double score = NaiveScore();
    if (score > bscore) {
        bscore = score;
        F0(i, n) {
            bans[i].first = pii(X1[i], Y1[i]);
            bans[i].second = pii(X2[i], Y2[i]);
        }
    }
}

void LoadBest() {
    F0(i, n) {
        X1[i] = bans[i].first.first;
        Y1[i] = bans[i].first.second;
        X2[i] = bans[i].second.first;
        Y2[i] = bans[i].second.second;
    }
    score = bscore;
}

int maxs = 1000;
const int d1 = 9;
const int d0 = d1+1;
const int d2 = 2*d0+1;
bool ExpandStep(int i) {
    while (X1[i] > bx1+d1 && d1+X2[i] < bx2 && Y1[i] > by1+d1 && d1+Y2[i] < by2 && (Y2[i] - Y1[i] + d2) * (X2[i] - X1[i] + d2) <= r[i]) {
        X1[i]-=d0;
        X2[i]+=d0;
        Y1[i]-=d0;
        Y2[i]+=d0;
    }

    int xlimit = min(r[i] / (Y2[i] - Y1[i] + 1) - (X2[i] - X1[i] + 1), maxs);
    int ylimit = min(r[i] / (X2[i] - X1[i] + 1) - (Y2[i] - Y1[i] + 1), maxs);
    if (xlimit == 0 && ylimit == 0) return false;

    biasind[0] = min(X1[i] - bx1, xlimit);
    biasind[1] = min(Y1[i] - by1, ylimit);
    biasind[2] = min(bx2 - X2[i], xlimit);
    biasind[3] = min(by2 - Y2[i], ylimit);

    F0(i, 4) biases[i] = (biasind[i] == 0) ? 0 : 1;
    int dir = biased(4);

    if (dir == -1) {
        return false;
    }
    switch (dir) {
    case 0: X1[i] -= biasind[0]; return true;
    case 1: Y1[i] -= biasind[1]; return true;
    case 2: X2[i] += biasind[2]; return true;
    case 3: Y2[i] += biasind[3]; return true;
    }
    return true;
}

void Expand(int i) {
    if ((bx2 - bx1 + 1) * (by2 - by1 + 1) <= r[i]) {
        X1[i] = bx1;
        Y1[i] = by1;
        X2[i] = bx2;
        Y2[i] = by2;
        return;
    }

    while (bx1 < X1[i] || X2[i] < bx2 || by1 < Y1[i] || Y2[i] < by2) {
        bool f = ExpandStep(i);
        if (!f) break;
    }
}

void SolveGreedy(bool scratch) {
    F0(i, M/CR+1) F0(j, M/CR+1) tn[i][j] = 0;
    F0(i, n) {
        coverx[i].n = 0;
        covery[i].n = 0;
    }
    if (!scratch) LoadBest(); else {
        F0(i, n) {
            X1[i] = X2[i] = X[i];
            Y1[i] = Y2[i] = Y[i];
        }
    }
    F0(i, n) AddT(i);
    F0(i, n) F0(j, n) if (i != j) {
        if (X1[i] <= X[j] && X2[i] >= X[j]) coverx[j].insert(i);
        if (Y1[i] <= Y[j] && Y2[i] >= Y[j]) covery[j].insert(i);
    }

    F0(i, n) id[i] = i;
    sort(id, id + n, [&](int i1, int i2) { return r[i1] < r[i2]; });

    vector<int> faker(n, 0);
    F0(i, n) faker[i] = r[i];
    int prep = 900;
    int itc = scratch ? prep : 700;
    F1(it, itc) {
        double x = 1.0 * it / prep;
        double p = pow(x, 0.5);
        F0(i, n) {
            if (it >= prep || !scratch) r[i] = faker[i];
            else r[i] = max(1, (int)(faker[i] * p));
            s[i] = (X2[i] - X1[i] + 1) * (Y2[i] - Y1[i] + 1);
        }
        sort(id, id + n, [&](int i1, int i2) { return 1LL*s[i1]*r[i2] < 1LL*s[i2]*r[i1]; });
        F0(u, n) {
            int i = id[u];
            qox1 = X1[i]; qox2 = X2[i]; qoy1 = Y1[i]; qoy2 = Y2[i];
            RemT(i);
            GenBoundaryFast(i);

            int changed = 0;
            Expand(i);
            ox1 = X1[i]; ox2 = X2[i]; oy1 = Y1[i]; oy2 = Y2[i]; oa = (X2[i]-X1[i]+1) * (Y2[i]-Y1[i]+1);
            if ((bx2-bx1+1) * (by2-by1+1) < r[i]) {
                X1[i] = X2[i] = X[i];
                Y1[i] = Y2[i] = Y[i];
                changed = 1;
                GenBoundary2(i, 4);
                Expand(i);
            }
            if (!changed) {
                int dir = 0;
                biasind[0] = min(X1[i] - bx1, X2[i] - X[i]);
                biasind[1] = min(Y1[i] - by1, Y2[i] - Y[i]);
                biasind[2] = min(bx2 - X2[i], X[i] - X1[i]);
                biasind[3] = min(by2 - Y2[i], Y[i] - Y1[i]);
                F0(i, 4) biases[i] = (biasind[i] == 0) ? 0 : 1;
                dir = biased(4);
                if (dir != -1 && my.nextInt(4)) {
                    if (dir == 0) X1[i] -= biasind[0], X2[i] -= biasind[0];
                    if (dir == 1) Y1[i] -= biasind[1], Y2[i] -= biasind[1];
                    if (dir == 2) X1[i] += biasind[2], X2[i] += biasind[2];
                    if (dir == 3) Y1[i] += biasind[3], Y2[i] += biasind[3];
                } else {
                    X1[i] = X2[i] = X[i];
                    Y1[i] = Y2[i] = Y[i];
                    GenBoundaryFaster(i);
                    Expand(i);
                }
            }
            double ll = 0.88;
            if ((X2[i]-X1[i]+1) * (Y2[i]-Y1[i]+1) < ll * oa) {
                X1[i] = ox1;
                Y1[i] = oy1;
                X2[i] = ox2;
                Y2[i] = oy2;
            }
            AddT(i);
            ChangeCovers(i);
        }
        if (it >= prep || !scratch) UpdateBest();
    }
    UpdateBest();
}

int bd[N][4];
int mv_out[N], mv_in[N];
set<pii> bar[N][4];
multiset<int> hb[N][4];
int st[N], sn;

bool DFS(int i) {
    if (sn > 20) return false;
    st[sn++] = i;
    Rect cur = nrect[i];

    int k = my.nextInt(4);
    int done = 0;
    F0(z, 4) {
        k = (k + 1) % 4;
        int hd = *hb[i][k].begin();
        if (!hd) continue;
        int room = RoomK(cur, k, r[i]);
        if (!room) continue;
        int ad = min(room, hd);

        int d = bar[i][k].begin()->first;

        // free growth
        // TODO: maybe still continue
        if (d != 0) {
            ad = min(ad, d);
            double new_delta = delta - upper(Area(cur), r[i]);
            int cd = ad;
            nrect[i] = Grow(cur, k, cd);
            new_delta += upper(Area(nrect[i]), r[i]);
            if (new_delta >= maxd) {
                delta = new_delta;
                return true;
            }
        } else if (!done) {
            // TODO: maybe we can do both, or maybe more
            auto it = bar[i][k].begin(); it++;
            ad = min(ad, it->first);
            if (ad == 0) continue;

            int j = bar[i][k].begin()->second;
            if (mv_out[j] != -1) continue;

            int cd = my.nextInt(1, ad);
            //cd = my.nextInt(cd, ad);
            //cd = my.nextInt(cd, ad);
            //cd = ad;

            mv_out[i] = k;

            double new_delta = delta - upper(Area(cur), r[i]) - upper(Area(rect[j]), r[j]);
            nrect[i] = Grow(cur, k, cd);
            nrect[j] = Grow(rect[j], opp(k), -cd);
            new_delta += upper(Area(nrect[i]), r[i]) + upper(Area(nrect[j]), r[j]);
            if (new_delta >= maxd) {
                delta = new_delta;
                st[sn++] = j;
                return true;
            }
            double old_delta = delta;
            delta = new_delta;
            done++;
            if (DFS(j)) {
                return true;
            }
            delta = old_delta;
            mv_out[i] = -1;
        }
    }

    sn--;
    return false;
}

void ToRect() {
    F0(i, n) rect[i] = Rect(pii(X1[i], Y1[i]), pii(X2[i], Y2[i]));
}

void FromRect() {
    F0(i, n) {
        X1[i] = rect[i].first.first;
        Y1[i] = rect[i].first.second;
        X2[i] = rect[i].second.first;
        Y2[i] = rect[i].second.second;
    }
}

void AddBarX(int i, int j) {
    int dist = rect[j].first.second - rect[i].second.second - 1;
    hb[i][3].insert(Y[j] - rect[i].second.second - 1);
    bar[i][3].insert(pii(dist, j));
    hb[j][1].insert(rect[j].first.second - Y[i] - 1);
    bar[j][1].insert(pii(dist, i));
}

void AddBarY(int i, int j) {
    int dist = rect[j].first.first - rect[i].second.first - 1;
    hb[i][2].insert(X[j] - rect[i].second.first - 1);
    bar[i][2].insert(pii(dist, j));
    hb[j][0].insert(rect[j].first.first - X[i] - 1);
    bar[j][0].insert(pii(dist, i));
}

void RemoveBarX(int i, int j) {
    int dist = rect[j].first.second - rect[i].second.second - 1;
    hb[i][3].erase(hb[i][3].find(Y[j] - rect[i].second.second - 1));
    bar[i][3].erase(pii(dist, j));
    hb[j][1].erase(hb[j][1].find(rect[j].first.second - Y[i] - 1));
    bar[j][1].erase(pii(dist, i));
}

void RemoveBarY(int i, int j) {
    int dist = rect[j].first.first - rect[i].second.first - 1;
    hb[i][2].erase(hb[i][2].find(X[j] - rect[i].second.first - 1));
    bar[i][2].erase(pii(dist, j));
    hb[j][0].erase(hb[j][0].find(rect[j].first.first - X[i] - 1));
    bar[j][0].erase(pii(dist, i));
}

void CalcBar() {
    F0(i, n) F0(k, 4) bar[i][k].clear();
    F0(i, n) {
        bar[i][0].insert(pii(rect[i].first.first, -1));
        bar[i][1].insert(pii(rect[i].first.second, -1));
        bar[i][2].insert(pii(M-1-rect[i].second.first, -1));
        bar[i][3].insert(pii(M-1-rect[i].second.second, -1));
        F0(k, 4) {
            hb[i][k].clear();
            hb[i][k].insert(bar[i][k].begin()->first);
        }
    }
    F0(i, n) F0(j, n) if (i != j) {
        if (IntersectsX(rect[i], rect[j]) && rect[i].first.second < rect[j].first.second) {
            AddBarX(i, j);
        }
        if (IntersectsY(rect[i], rect[j]) && rect[i].first.first < rect[j].first.first) {
            AddBarY(i, j);
        }
    }
}

void RemoveBar(int i) {
    F0(j, n) if (lck[j]) {
        if (IntersectsX(rect[i], rect[j])) {
            if (rect[i].first.second < rect[j].first.second) RemoveBarX(i, j); else RemoveBarX(j, i);
        }
        if (IntersectsY(rect[i], rect[j])) {
            if (rect[i].first.first < rect[j].first.first) RemoveBarY(i, j); else RemoveBarY(j, i);
        }
    }
    F0(k, 4) {
        hb[i][k].clear();
        bar[i][k].clear();
    }
}

void AddBar(int i) {
    bar[i][0].insert(pii(rect[i].first.first, -1));
    bar[i][1].insert(pii(rect[i].first.second, -1));
    bar[i][2].insert(pii(M-1-rect[i].second.first, -1));
    bar[i][3].insert(pii(M-1-rect[i].second.second, -1));
    F0(k, 4) {
        hb[i][k].insert(bar[i][k].begin()->first);
    }
    F0(j, n) if (lck[j]) {
        if (IntersectsX(rect[i], rect[j])) {
            if (rect[i].first.second < rect[j].first.second) AddBarX(i, j); else AddBarX(j, i);
        }
        if (IntersectsY(rect[i], rect[j])) {
            if (rect[i].first.first < rect[j].first.first) AddBarY(i, j); else AddBarY(j, i);
        }
    }
}

void FindPaths(int i) {
    F0(j, n) mv_out[j] = -1;
    nrect[i] = rect[i];
    delta = 0.0;
    if (DFS(i)) {
        if (sn == 0) throw;
        return;
    }
    if (sn != 0) throw;
}

void SolveSA() {
    PR(Elapsed());
    score = 0.0;
    F0(i, n) {
        s[i] = (X2[i] - X1[i] + 1) * (Y2[i] - Y1[i] + 1);
        score += upper(s[i], r[i]);
    }
    int numStages = 20;
    //double T1 = 0.0, T2 = 0.1000000, T = 0.0;
    double T1 = 0.0, T2 = 0.000010, T = 0.0;
    int itc = 1000000000;
    int tot = 1, acc = 0;
    int lastStage = -1;

    ToRect();
    CalcBar();

    double startTime = Elapsed();
    double endTime = TL2;

    if (startTime >= endTime) return;

    F0(iter, itc) {
        if ((iter & 255) == 0) {
            //double r = 1.0 * (itc - iter) / itc;
            double r = (endTime - Elapsed()) / (endTime - startTime);
            if (1) {
                int newStage = (1 - r) * numStages;
                if (newStage != lastStage) {
                    lastStage = newStage;
                    cerr << "St: " << newStage << " " << 100.0*acc/tot << " " << score << endl;
                }
            }
            if (r <= 0) {
                PR(iter);
                break;
            }
            T = T1 + (T2 - T1) * r;
        }

        int i = -1;
        maxd = 0.0;
        sn = 0;
        i = my.nextInt(n);
        if (MaxedOut(rect[i], r[i])) {
            //if (my.nextInt(8)) continue;
            int k = my.nextInt(4);
            int hd = *hb[i][k].begin();
            int d = bar[i][k].begin()->first;
            int ad = min(hd, d);
            int back = 0;
            switch (k) {
            case 0: back = rect[i].second.first - X[i]; break;
            case 1: back = rect[i].second.second - Y[i]; break;
            case 2: back = X[i] - rect[i].first.first; break;
            case 3: back = Y[i] - rect[i].first.second; break;
            };
            ad = min(ad, back);
            if (ad > 0) {
                int cd = my.nextInt(1, ad);
                nrect[i] = Grow(rect[i], k, cd);
                nrect[i] = Grow(nrect[i], opp(k), -cd);
                sn = 1;
                st[0] = i;
                delta = 0.0;
            }
        } else {
            maxd = -T * logs[my.nextInt(LOGN)];
            tot++;
            FindPaths(i);
        }

        if (sn == 0) continue;

        bool unl = false;
        F0(j, sn) {
            int u = st[j];
            F0(k, j) {
                int v = st[k];
                if (Intersects(nrect[u], nrect[v])) {
                    unl = true;
                    j = sn;
                }
            }
            if (!Contains(nrect[u], X[u], Y[u])) {
                throw;
            }
        }
        if (unl) continue;

        if (delta < 0.0) {
            acc++;
        }
        score += delta;

        F0(i, n) lck[i] = 1;
        F0(j, sn) {
            int u = st[j];
            lck[u] = 0;
            RemoveBar(u);
        }
        F0(j, sn) {
            int u = st[j];
            rect[u] = nrect[u];
        }
        F0(j, sn) {
            int u = st[j];
            AddBar(u);
            lck[u] = 1;
        }

        if (score > bscore) {
            FromRect();
            UpdateBest();
        }

        FromRect();
    }
}

vector<int> active;
void SolveLocal() {
    F0(i, M/CR+1) F0(j, M/CR+1) tn[i][j] = 0;
    F0(i, n) {
        coverx[i].n = 0;
        covery[i].n = 0;
    }
    F0(i, n) AddT(i);
    F0(i, n) F0(j, n) if (i != j) {
        if (X1[i] <= X[j] && X2[i] >= X[j]) coverx[j].insert(i);
        if (Y1[i] <= Y[j] && Y2[i] >= Y[j]) covery[j].insert(i);
    }

    int nn = SZ(active);
    F0(i, nn) id[i] = active[i];

    double oscore = 0.0;
    F0(i, n) {
        s[i] = (X2[i] - X1[i] + 1) * (Y2[i] - Y1[i] + 1);
        oscore += upper(s[i], r[i]);
    }
    F0(u, nn) {
        int i = id[u];
        s[i] = (X2[i] - X1[i] + 1) * (Y2[i] - Y1[i] + 1);
        oscore -= upper(s[i], r[i]);
    }

    vector<int> faker(n, 0);
    F0(i, n) faker[i] = r[i];
    int prep = 50;
    int itc = 75;
    F1(it, itc) {
        double x = 1.0 * it / prep;
        double p = pow(x, 0.5);
        score = oscore;
        F0(u, nn) {
            int i = id[u];
            if (it >= prep) r[i] = faker[i];
            else r[i] = max(1, (int)(faker[i] * p));
            s[i] = (X2[i] - X1[i] + 1) * (Y2[i] - Y1[i] + 1);
            if (it >= prep) score += upper(s[i], r[i]);
        }
        sort(id, id + nn, [&](int i1, int i2) { return 1LL*s[i1]*r[i2] < 1LL*s[i2]*r[i1]; });
        F0(u, nn) {
            int i = id[u];
            if (s[i] > r[i]) continue;
            if (it >= prep) score -= upper(s[i], r[i]);
            qox1 = X1[i]; qox2 = X2[i]; qoy1 = Y1[i]; qoy2 = Y2[i];
            RemT(i);
            GenBoundaryFast(i);

            int changed = 0;
            Expand(i);
            ox1 = X1[i]; ox2 = X2[i]; oy1 = Y1[i]; oy2 = Y2[i]; oa = (X2[i]-X1[i]+1) * (Y2[i]-Y1[i]+1);
            if ((bx2-bx1+1) * (by2-by1+1) < r[i]) {
                X1[i] = X2[i] = X[i];
                Y1[i] = Y2[i] = Y[i];
                changed = 1;
                GenBoundary2(i, 2);
                //GenBoundaryFaster(i);
                Expand(i);
            }
            if (!changed) {
                int dir = 0;
                biasind[0] = min(X1[i] - bx1, X2[i] - X[i]);
                biasind[1] = min(Y1[i] - by1, Y2[i] - Y[i]);
                biasind[2] = min(bx2 - X2[i], X[i] - X1[i]);
                biasind[3] = min(by2 - Y2[i], Y[i] - Y1[i]);
                F0(i, 4) biases[i] = (biasind[i] == 0) ? 0 : 1;
                dir = biased(4);
                if (dir != -1 && my.nextInt(4)) {
                    if (dir == 0) X1[i] -= biasind[0], X2[i] -= biasind[0];
                    if (dir == 1) Y1[i] -= biasind[1], Y2[i] -= biasind[1];
                    if (dir == 2) X1[i] += biasind[2], X2[i] += biasind[2];
                    if (dir == 3) Y1[i] += biasind[3], Y2[i] += biasind[3];
                } else {
                    X1[i] = X2[i] = X[i];
                    Y1[i] = Y2[i] = Y[i];
                    GenBoundaryFaster(i);
                    Expand(i);
                }
            }
            double ll = 0.88;
            if ((X2[i]-X1[i]+1) * (Y2[i]-Y1[i]+1) < ll * oa) {
                X1[i] = ox1;
                Y1[i] = oy1;
                X2[i] = ox2;
                Y2[i] = oy2;
            }
            AddT(i);
            ChangeCovers(i);

            if (it >= prep) {
                s[i] = (X2[i] - X1[i] + 1) * (Y2[i] - Y1[i] + 1);
                score += upper(s[i], r[i]);
                if (score > bscore) UpdateBest();
            }
            //if (it >= prep) UpdateBest();
        }
    }
}


vector<pii> adj[N];
void LocalS() {
    F0(i, n) {
        F0(j, n) if (j != i) adj[i].push_back(make_pair(abs(X[i]-X[j]) + abs(Y[i]-Y[j]), j));
        sort(adj[i].begin(), adj[i].end());
    }
    while (Elapsed() < TL1) {
        c_o[2]++;

        active.clear();
        int MINT = 11;
        int MAXT = 15;

        int T = my.nextInt(MINT, MAXT);
        int st = my.nextInt(n);

        set<int> chosen;
        chosen.insert(st);
        F0(j, T) {
            chosen.insert(adj[st][j].second);
        }

        LoadBest();
        if (1) {
            SUBX1 = M; SUBX2 = 0;
            SUBY1 = M; SUBY2 = 0;
            for (int i : chosen) {
                SUBX1 = min(SUBX1, X1[i] - 1000);
                SUBY1 = min(SUBY1, Y1[i] - 1000);
                SUBX2 = max(SUBX2, X2[i] + 1000);
                SUBY2 = max(SUBY2, Y2[i] + 1000);
            }
            SUBX1 = max(SUBX1, 0);
            SUBY1 = max(SUBY1, 0);
            SUBX2 = min(SUBX2, M-1);
            SUBY2 = min(SUBY2, M-1);
        }

        active.clear();
        for (int i : chosen) //if (my.nextInt(3))
            active.push_back(i);
        double pscore = bscore;
        for (int i : active) {
            X1[i] = X2[i] = X[i];
            Y1[i] = Y2[i] = Y[i];
        }
        SolveLocal();
        //cerr << pscore << " " << NaiveScore() << endl;
        if (bscore > pscore) {
            //PR(bscore - pscore);
            //return;
        }
    }
}

void Solve(bool print_ans = true) {
    Init();
    PR(n);
    F0(i, LOGN) logs[i] = -log((i+0.5)/LOGN);

    vector<pii> sorted;
    X[n] = inf; Y[n] = inf;
    F0(i, n + 1) sorted.push_back(pii(X[i], i));
    sort(sorted.begin(), sorted.end());
    F0(i, n) nx[sorted[i].second] = sorted[i + 1].second;
    int j = 0;
    F0(i, M) {
        while (i > sorted[j].first) j++;
        fx[i] = sorted[j].second;
    }
    F0(i, n + 1) sorted[i] = pii(Y[i], i);
    sort(sorted.begin(), sorted.end());
    F0(i, n) ny[sorted[i].second] = sorted[i + 1].second;
    j = 0;
    F0(i, M) {
        while (i > sorted[j].first) j++;
        fy[i] = sorted[j].second;
    }

    score = 0.0;
    F0(i, n) {
        X1[i] = X2[i] = X[i];
        Y1[i] = Y2[i] = Y[i];
        s[i] = 1;
        score += upper(s[i], r[i]);
    }
    UpdateBest();

    SUBX1 = 0; SUBX2 = M-1; SUBY1 = 0; SUBY2 = M-1;

    if (1) {
        int lc = 0;
        while (Elapsed() < TL0) {
            lc++;
            if (lc % 10 == 0) PR(lc);
            SolveGreedy(true);
        }
        PR(lc);
        SolveGreedy(false);
        LocalS();
        PR(c_o[1]);
        LoadBest();
    }

    if (1) {
        SolveSA();
        LoadBest();
    }

    F0(i, n) {
        if (print_ans) cout << X1[i] << " " << Y1[i] << " " << X2[i]+1 << " " << Y2[i]+1 << endl;
        if (X1[i] < 0 || X1[i] > X2[i] || X2[i] >= M || Y1[i] < 0 || Y1[i] > Y2[i] || Y2[i] >= M) {
            cerr << "Bad rectangles" << endl;
            throw;
        }
        F0(j, i) {
            if (X1[i] > X2[j] || X1[j] > X2[i] || Y1[i] > Y2[j] || Y1[j] > Y2[i]) continue;
            cerr << i << " and " << j << " intersect" << endl;
            throw;
        }
        double rw = 0;
        if (X1[i] <= X[i] && X[i] <= X2[i] && Y1[i] <= Y[i] && Y[i] <= Y2[i]) {
            rw = raw((X2[i]-X1[i]+1) * (Y2[i]-Y1[i]+1), r[i]);
        }
        act_score += rw;
    }    //if (!print_ans)
}












void RunVis() {
    cin >> n;
    F0(i, n) {
        cin >> X[i] >> Y[i] >> r[i];
    }
    Solve();
}

void RunSeed(int seed) {
    Rand rnd(seed);
    n = 50 * pow(4.0, rnd.nextDouble()) + 0.5;
    set<pii> S;
    F0(i, n) {
        while (1) {
            X[i] = rnd.nextInt(10000);
            Y[i] = rnd.nextInt(10000);
            if (!S.count(pii(X[i], Y[i]))) {
                S.insert(pii(X[i], Y[i]));
                break;
            }
        }
    }
    set<int> Q;
    Q.insert(0);
    Q.insert(M*M);
    F0(i, n-1) {
        while (1) {
            int q = rnd.nextInt(M*M-1) + 1;
            if (!Q.count(q)) {
                Q.insert(q);
                break;
            }
        }
    }
    auto it = Q.begin();
    F0(i, n) {
        auto it2 = it;
        it2++;
        r[i] = *it2 - *it;
        it = it2;
    }
    Solve(false);
}


int main(int argc, char* argv[]) {
#ifdef RED
    //auto x = freopen("log.txt", "w", stderr); (void)x;
#endif
    if (argc > 1) {
        int seed = atoi(argv[1]);
        if (seed == 0) {
            ignore = freopen("x.in", "r", stdin);
            ignore = freopen("x.out", "w", stdout);
            RunVis();
        } else {
            PR(seed);
            RunSeed(seed);
        }
    } else {
        //ignore = freopen("x.out", "w", stdout);
        RunVis();
    }

    return 0;
}
