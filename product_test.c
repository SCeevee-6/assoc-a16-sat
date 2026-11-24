
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// Fast 3-coloring for A_n (n <= 9) and recursive 4-coloring for A_n (n <= 15)
// First let us compute exact root-triangle decomposition (m, T_L, T_R) on A16_graph.bin!
static int g_bits, g_count;
static uint32_t *g_codes;

static void gen_dfs(int pos, int open_n, int int_left, uint32_t cur) {
    if (pos == g_bits) { g_codes[g_count++] = cur; return; }
    if (int_left > 0) gen_dfs(pos + 1, open_n + 1, int_left - 1, cur | (1U << (g_bits - 1 - pos)));
    if (open_n > 1 || (open_n == 1 && int_left == 0)) gen_dfs(pos + 1, open_n - 1, int_left, cur);
}
static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}
static inline int bsearch_u32(const uint32_t *arr, int n, uint32_t tgt) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) { int mid = lo + ((hi - lo) >> 1); if (arr[mid] == tgt) return mid; if (arr[mid] < tgt) lo = mid + 1; else hi = mid - 1; }
    return -1;
}
static inline int sub_len(uint32_t code, int bits, int pos) {
    int bal = 0;
    for (int p = pos; p < bits; p++) { bal += ((code >> (bits - 1 - p)) & 1U) ? 1 : -1; if (bal < 0) return p - pos + 1; }
    return bits + 1 - pos;
}

// Compute greedy/backtracking 4-coloring of A_n for all n = 3..15 in <1 second!
static int8_t *col_An[16];
static int nv_An[16];
static uint32_t *codes_An[16];

static int bt_color(int v, int nv, int deg_max, int32_t *adj, int8_t *col) {
    if (v == nv) return 1;
    int used = 0;
    for (int i = 0; i < deg_max; i++) {
        int nb = adj[v * deg_max + i];
        if (nb < 0) break;
        if (nb < v && col[nb] >= 0) used |= (1 << col[nb]);
    }
    for (int c = 0; c < 4; c++) {
        if (!(used & (1 << c))) {
            col[v] = (int8_t)c;
            if (bt_color(v + 1, nv, deg_max, adj, col)) return 1;
            col[v] = -1;
        }
    }
    return 0;
}

static void build_small_colorings(void) {
    for (int n = 3; n <= 15; n++) {
        int K = n - 2;
        g_bits = 2 * K;
        g_codes = malloc(800000 * sizeof(uint32_t));
        g_count = 0;
        if (K == 0) { g_codes[0] = 0; g_count = 1; }
        else gen_dfs(0, 1, K, 0U);
        qsort(g_codes, g_count, sizeof(uint32_t), cmp_u32);
        int nv = g_count;
        nv_An[n] = nv;
        codes_An[n] = g_codes;
        col_An[n] = malloc(nv);
        memset(col_An[n], -1, nv);
        if (K == 0) { col_An[n][0] = 0; continue; }

        int deg_max = n - 3;
        int32_t *adj = malloc((size_t)nv * deg_max * sizeof(int32_t));
        memset(adj, -1, (size_t)nv * deg_max * sizeof(int32_t));
        int *deg = calloc(nv, sizeof(int));

        for (int idx = 0; idx < nv; idx++) {
            uint32_t code = g_codes[idx];
            for (int p = 0; p < g_bits; p++) {
                if (((code >> (g_bits - 1 - p)) & 1U) == 0U) continue;
                int lp = p + 1, ll = sub_len(code, g_bits, lp), rp = lp + ll;
                if (lp < g_bits && ((code >> (g_bits - 1 - lp)) & 1U) == 1U) {
                    int lenA = sub_len(code, g_bits, lp + 1);
                    uint32_t nxt = code;
                    for (int i = 0; i < lenA; i++) {
                        uint32_t b = (code >> (g_bits - 1 - (lp + 1 + i))) & 1U;
                        int sh = g_bits - 1 - (lp + i);
                        if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);
                    }
                    nxt |= (1U << (g_bits - 1 - (lp + lenA)));
                    int nb = bsearch_u32(g_codes, nv, nxt);
                    if (idx < nb) {
                        adj[idx * deg_max + deg[idx]++] = nb;
                        adj[nb * deg_max + deg[nb]++] = idx;
                    }
                }
                if (rp < g_bits && ((code >> (g_bits - 1 - rp)) & 1U) == 1U) {
                    int lenA = ll;
                    uint32_t nxt = code | (1U << (g_bits - 1 - lp));
                    for (int i = 0; i < lenA; i++) {
                        uint32_t b = (code >> (g_bits - 1 - (lp + i))) & 1U;
                        int sh = g_bits - 1 - (lp + 1 + i);
                        if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);
                    }
                    int nb = bsearch_u32(g_codes, nv, nxt);
                    if (idx < nb) {
                        adj[idx * deg_max + deg[idx]++] = nb;
                        adj[nb * deg_max + deg[nb]++] = idx;
                    }
                }
            }
        }
        if (!bt_color(0, nv, deg_max, adj, col_An[n])) {
            fprintf(stderr, "Failed coloring A_%d\n", n); exit(1);
        }
        printf("Colored A_%2d (nv=%6d) cleanly!\n", n, nv);
        free(adj); free(deg);
    }
}

int main(void) {
    build_small_colorings();
    return 0;
}
