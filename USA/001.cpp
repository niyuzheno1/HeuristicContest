const double TL2 = 9.7;
const double TL1 = 4*TL2/6.0;
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
Rand my(2020);
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
double bscore, lscore;
Rect bans[N], lans[N];
double delta, act_score, score;
int bad, dir, banned;
int X1[N], Y1[N], X2[N], Y2[N], s[N];
int SUBX1, SUBX2, SUBY1, SUBY2, SPLITX;
ll sumbounded;
const int LOGN = 1 << 16;
double logs[LOGN];

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

const int CR = 1024;
const int MAXG = M/CR+10;
int tn[MAXG][MAXG], t[MAXG][MAXG][128];
int BX1, BX2, BY1, BY2;
int px1, px2, py1, py2;
int bx1, bx2, by1, by2;
int ox1, oy1, ox2, oy2, oa;
int qox1, qoy1, qox2, qoy2;
int fx[M], fy[M], nx[N], ny[N];
int undecided[N], unn;

void PrepG(int i) {
    BX1 = X1[i] >> 10;
    BY1 = Y1[i] >> 10;
    BX2 = X2[i] >> 10;
    BY2 = Y2[i] >> 10;
}

void PrepG() {
    BX1 = bx1 >> 10;
    BY1 = by1 >> 10;
    BX2 = bx2 >> 10;
    BY2 = by2 >> 10;
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
    throw;
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
    ox1 = qox1;
    ox2 = qox2;
    oy1 = qoy1;
    oy2 = qoy2;
    if (X1[i] != ox1) {
        if (X1[i] < ox1) {
            for (int j = fx[X1[i]]; X[j] < ox1; j = nx[j]) coverx[j].insert(i);
        } else {
            for (int j = fx[ox1]; X[j] < X1[i]; j = nx[j]) coverx[j].erase(i);
        }
    }
    if (X2[i] != ox2) {
        if (X2[i] > ox2) {
            for (int j = fx[ox2+1]; X[j] <= X2[i]; j = nx[j]) coverx[j].insert(i);
        } else {
            for (int j = fx[X2[i]+1]; X[j] <= ox2; j = nx[j]) coverx[j].erase(i);
        }
    }
    if (Y1[i] != oy1) {
        if (Y1[i] < oy1) {
            for (int j = fy[Y1[i]]; Y[j] < oy1; j = ny[j]) covery[j].insert(i);
        } else {
            for (int j = fy[oy1]; Y[j] < Y1[i]; j = ny[j]) covery[j].erase(i);
        }
    }
    if (Y2[i] != oy2) {
        if (Y2[i] > oy2) {
            for (int j = fy[oy2+1]; Y[j] <= Y2[i]; j = ny[j]) covery[j].insert(i);
        } else {
            for (int j = fy[Y2[i]+1]; Y[j] <= oy2; j = ny[j]) covery[j].erase(i);
        }
    }
}

int seenZ[N], Z;
void GenBoundaryFast(int i) {
    AT(9);
    c_o[9]++;
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
                c_o[2]++;
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
        c_o[3]++;
        int j = undecided[uuu];
        if (X1[j] > X2[i] && (dx = bx2-X1[j]+1) > 0) {
            if (Y1[j] > Y2[i] && (dy = by2-Y1[j]+1) > 0) {
                c_o[2]++;
                if (my.nextInt(dx+dy) < dx) {
                    by2 -= dy;
                } else {
                    bx2 -= dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                c_o[2]++;
                if (my.nextInt(dx+dy) < dx) {
                    by1 += dy;
                } else {
                    bx2 -= dx;
                }
            }
        }
        if (X1[j] < X2[i] && (dx = X2[j]+1-bx1)>0) {
            if (Y1[j] > Y2[i] && (dy = by2-Y1[j]+1) > 0) {
                c_o[2]++;
                if (my.nextInt(dx+dy) < dx) {
                    by2 -= dy;
                } else {
                    bx1 += dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                c_o[2]++;
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
    AT(8);
    c_o[8]++;
    bx1 = SUBX1; by1 = SUBY1;
    bx2 = SUBX2; by2 = SUBY2;

    unn = 0;
    Z++;
    F0(u, coverx[i].n) {
        c_o[7]++;
        int j = coverx[i].a[u];
        seenZ[j] = Z;
        if (Y2[j] < Y[i]) by1 = max(by1, Y2[j]+1);
        else by2 = min(by2, Y1[j]-1);
    }
    F0(u, covery[i].n) {
        int j = covery[i].a[u];
        c_o[7]++;
        seenZ[j] = Z;
        if (X2[j] < X[i]) bx1 = max(bx1, X2[j]+1);
        else bx2 = min(bx2, X1[j]-1);
    }

    PrepG();
    for (int x = BX1; x <= BX2; x++) for (int y = BY1; y <= BY2; y++) {
        F0(k, tn[x][y]) {
            int j = t[x][y][k];
            if (seenZ[j] != Z) {
                c_o[2]++;
                seenZ[j] = Z;
                undecided[unn++] = j;
            }
        }
    }

    int dx = 0, dy = 0;
    F0(uuu, unn) {
        c_o[3]++;
        int j = undecided[uuu];
        if (X1[j] > X2[i] && (dx = bx2-X1[j]+1) > 0) {
            if (Y1[j] > Y2[i] && (dy = by2-Y1[j]+1) > 0) {
                c_o[2]++;
                biases[0] = dx * (by2 - by1 + 1);
                biases[1] = dy * (bx2 - bx1 + 1);
                if (my.nextInt(biases[0]+biases[1]) < biases[0]) {
                    by2 -= dy;
                } else {
                    bx2 -= dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                c_o[2]++;
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
                c_o[2]++;
                biases[0] = dx * (by2 - by1 + 1);
                biases[1] = dy * (bx2 - bx1 + 1);
                if (my.nextInt(biases[0]+biases[1]) < biases[0]) {
                    by2 -= dy;
                } else {
                    bx1 += dx;
                }
            }
            if (Y1[j] < Y2[i] && (dy = Y2[j]+1-by1) > 0) {
                c_o[2]++;
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
    AT(1);
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
    lscore = 0.0;
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

    SUBX1 = 0; SUBX2 = M-1; SUBY1 = 0; SUBY2 = M-1;

    vector<int> faker(n, 0);
    F0(i, n) faker[i] = r[i];
    int prep = 1000;
    F0(it, 1000) {
        F0(i, n) {
            if (it >= prep || !scratch) r[i] = faker[i];
            else {
                double x = 1.0 * (it + 1) / prep;
                double p = pow(x, 0.5);
                r[i] = max(1, (int)(1LL * faker[i] * p));
            }
        }
        F0(i, n) s[i] = (X2[i] - X1[i] + 1) * (Y2[i] - Y1[i] + 1);
        sort(id, id + n, [&](int i1, int i2) { return 1.0*s[i1]/r[i1] < 1.0*s[i2]/r[i2]; });
        F0(u, n) {
            int i = id[u];

            qox1 = X1[i]; qox2 = X2[i]; qoy1 = Y1[i]; qoy2 = Y2[i];

            RemT(i);

            GenBoundaryFast(i);

            int changed = 0;
            Expand(i);
            ox1 = X1[i]; ox2 = X2[i]; oy1 = Y1[i]; oy2 = Y2[i]; oa = (X2[i]-X1[i]+1) * (Y2[i]-Y1[i]+1);
            if ((bx2-bx1+1) * (by2-by1+1) <= r[i]) {
                X1[i] = X2[i] = X[i];
                Y1[i] = Y2[i] = Y[i];
                changed = 1;
                GenBoundaryFaster(i);
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
                ChangeCovers(i);
            } else {
                ChangeCovers(i);
            }
            AddT(i);
        }
        if (it >= prep || !scratch) UpdateBest();
    }
    UpdateBest();
}


void Solve() {
    Init();
    PR(n);
    F0(i, LOGN) logs[i] = -log((i+0.5)/LOGN);

    if (0) {
        int maxr = 0;
        F0(i, n) {
            maxr = max(maxr, r[i]);
        }
        PR(maxr);
    }

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

    if (1) {
        int itc = 16;
        F0(u, itc) {
            SolveGreedy(u < itc - 1);
        }
        LoadBest();
    }

    F0(i, n) {
        cout << X1[i] << " " << Y1[i] << " " << X2[i]+1 << " " << Y2[i]+1 << endl;
    }
    Report();
}













void RunVis() {
    cin >> n;
    F0(i, n) {
        cin >> X[i] >> Y[i] >> r[i];
    }
    Solve();
}


int main() {
    RunVis();

    return 0;
}
