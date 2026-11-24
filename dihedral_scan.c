
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void code_to_diags(uint32_t code, int n_gon, int diags[16][2]) {
    int node_type[64], node_l[64], node_r[64];
    int nsp = 0, leaf = 0, d_idx = 0, bits = 2 * (n_gon - 2);
    for (int p = 0; p < bits; p++) {
        int bit = (int)((code >> (bits - 1 - p)) & 1U);
        node_type[nsp] = bit;
        if (bit == 0) { node_l[nsp] = leaf; node_r[nsp] = leaf + 1; leaf++; }
        nsp++;
        while (nsp >= 3 && node_type[nsp - 3] == 1 && node_type[nsp - 2] == 0 && node_type[nsp - 1] == 0) {
            int l = node_l[nsp - 2], r = node_r[nsp - 1];
            if (!(l == 0 && r == n_gon - 1)) { diags[d_idx][0] = l; diags[d_idx][1] = r; d_idx++; }
            nsp -= 2; node_type[nsp - 1] = 0; node_l[nsp - 1] = l; node_r[nsp - 1] = r;
        }
    }
    node_type[nsp] = 0; node_l[nsp] = leaf; node_r[nsp] = leaf + 1; nsp++;
    while (nsp >= 3 && node_type[nsp - 3] == 1 && node_type[nsp - 2] == 0 && node_type[nsp - 1] == 0) {
        int l = node_l[nsp - 2], r = node_r[nsp - 1];
        if (!(l == 0 && r == n_gon - 1)) { diags[d_idx][0] = l; diags[d_idx][1] = r; d_idx++; }
        nsp -= 2; node_type[nsp - 1] = 0; node_l[nsp - 1] = l; node_r[nsp - 1] = r;
    }
}

static uint32_t diags_to_code(int n_gon, int diags[16][2]) {
    uint8_t adj[16][16]; memset(adj, 0, sizeof(adj));
    for (int i = 0; i < n_gon - 3; i++) { int a = diags[i][0], b = diags[i][1]; adj[a][b] = adj[b][a] = 1; }
    uint32_t code = 0; int pos = 0, bits = 2 * (n_gon - 2), sl[64], sr[64], sp = 0;
    sl[sp] = 0; sr[sp] = n_gon - 1; sp++;
    while (sp > 0) {
        sp--; int l = sl[sp], r = sr[sp];
        if (r - l == 1) { if (pos < bits) pos++; }
        else {
            int mid = -1;
            for (int k = l + 1; k < r; k++) if (((k - l == 1) || adj[l][k]) && ((r - k == 1) || adj[k][r])) { mid = k; break; }
            code |= (1U << (bits - 1 - pos)); pos++;
            sl[sp] = mid; sr[sp] = r; sp++; sl[sp] = l; sr[sp] = mid; sp++;
        }
    }
    return code;
}

static int bsearch_u32(const uint32_t *arr, int n, uint32_t tgt) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) { int mid = lo + ((hi - lo) >> 1); if (arr[mid] == tgt) return mid; if (arr[mid] < tgt) lo = mid + 1; else hi = mid - 1; }
    return -1;
}

int main(void) {
    FILE *fg = fopen("A16_graph.bin", "rb");
    int32_t nv; int64_t ne;
    fread(&nv, sizeof(int32_t), 1, fg); fread(&ne, sizeof(int64_t), 1, fg);
    uint32_t *codes = malloc(nv * sizeof(uint32_t));
    int32_t *eu = malloc(ne * sizeof(int32_t)), *ev = malloc(ne * sizeof(int32_t));
    fread(codes, sizeof(uint32_t), nv, fg); fread(eu, sizeof(int32_t), ne, fg); fread(ev, sizeof(int32_t), ne, fg);
    fclose(fg);

    int32_t *rot1 = malloc(nv * sizeof(int32_t));
    int32_t *ref0 = malloc(nv * sizeof(int32_t));
    int32_t *ref1 = malloc(nv * sizeof(int32_t));
    for (int i = 0; i < nv; i++) {
        int d0[16][2], d1[16][2], d2[16][2];
        code_to_diags(codes[i], 16, d0);
        for (int d = 0; d < 13; d++) {
            int a = (d0[d][0] + 1) & 15, b = (d0[d][1] + 1) & 15;
            if (a > b) { int t = a; a = b; b = t; } d1[d][0] = a; d1[d][1] = b;
            a = (15 - d0[d][0]) & 15; b = (15 - d0[d][1]) & 15;
            if (a > b) { int t = a; a = b; b = t; } d2[d][0] = a; d2[d][1] = b;
        }
        rot1[i] = bsearch_u32(codes, nv, diags_to_code(16, d1));
        ref0[i] = bsearch_u32(codes, nv, diags_to_code(16, d2));
        for (int d = 0; d < 13; d++) {
            int a = (16 - d0[d][0]) & 15, b = (16 - d0[d][1]) & 15;
            if (a > b) { int t = a; a = b; b = t; } d2[d][0] = a; d2[d][1] = b;
        }
        ref1[i] = bsearch_u32(codes, nv, diags_to_code(16, d2));
    }

    int S4[24][4]; int sp = 0;
    for (int a=0;a<4;a++) for (int b=0;b<4;b++) if (b!=a) for (int c=0;c<4;c++) if (c!=a&&c!=b) for (int d=0;d<4;d++) if (d!=a&&d!=b&&d!=c) {
        S4[sp][0]=a; S4[sp][1]=b; S4[sp][2]=c; S4[sp][3]=d; sp++;
    }

    int ks[] = {8, 4, 2};
    int32_t *orb = malloc(nv * sizeof(int32_t));
    uint8_t *col_map = malloc(nv);
    int q[64];

    for (int ki = 0; ki < 3; ki++) {
        int k = ks[ki], step = 16 / k;
        for (int rtype = 0; rtype < 2; rtype++) {
            int32_t *ref = (rtype == 0) ? ref0 : ref1;
            for (int pr = 0; pr < 24; pr++) {
                int ok_r = 1;
                for (int c = 0; c < 4; c++) {
                    int cur = c; for (int s = 0; s < k; s++) cur = S4[pr][cur];
                    if (cur != c) ok_r = 0;
                }
                if (!ok_r) continue;
                for (int pf = 0; pf < 24; pf++) {
                    int ok_f = 1;
                    for (int c = 0; c < 4; c++) {
                        if (S4[pf][S4[pf][c]] != c) ok_f = 0;
                        int lhs = S4[pf][S4[pr][S4[pf][c]]];
                        if (S4[pr][lhs] != c) ok_f = 0;
                    }
                    if (!ok_f) continue;

                    memset(orb, -1, nv * sizeof(int32_t));
                    int norb = 0, consistent = 1;
                    for (int i = 0; i < nv && consistent; i++) {
                        if (orb[i] != -1) continue;
                        orb[i] = norb; col_map[i] = 0;
                        int qh = 0, qt = 0; q[qt++] = i;
                        while (qh < qt && consistent) {
                            int u = q[qh++], p_u = col_map[u];
                            int v1 = u; for (int s = 0; s < step; s++) v1 = rot1[v1];
                            int p_v1 = -1;
                            for (int idx = 0; idx < 24; idx++) {
                                if (S4[idx][0]==S4[pr][S4[p_u][0]] && S4[idx][1]==S4[pr][S4[p_u][1]] && S4[idx][2]==S4[pr][S4[p_u][2]]) { p_v1 = idx; break; }
                            }
                            if (orb[v1] == -1) { orb[v1] = norb; col_map[v1] = (uint8_t)p_v1; q[qt++] = v1; }
                            else if (col_map[v1] != p_v1) consistent = 0;

                            int v2 = ref[u];
                            int p_v2 = -1;
                            for (int idx = 0; idx < 24; idx++) {
                                if (S4[idx][0]==S4[pf][S4[p_u][0]] && S4[idx][1]==S4[pf][S4[p_u][1]] && S4[idx][2]==S4[pf][S4[p_u][2]]) { p_v2 = idx; break; }
                            }
                            if (orb[v2] == -1) { orb[v2] = norb; col_map[v2] = (uint8_t)p_v2; q[qt++] = v2; }
                            else if (col_map[v2] != p_v2) consistent = 0;
                        }
                        norb++;
                    }
                    if (!consistent) continue;

                    int self_loop = 0;
                    for (int64_t e = 0; e < ne; e++) {
                        int u = eu[e], v = ev[e];
                        if (orb[u] == orb[v] && col_map[u] == col_map[v]) { self_loop = 1; break; }
                    }
                    if (!self_loop) {
                        printf("SURVIVED: D_%d (ref%d) | rot=(%d%d%d%d) ref=(%d%d%d%d) -> %d orbits!\n",
                               S4[pf][0], S4[pf][1], S4[pf][2], S4[pf][3], norb);
                    }
                }
            }
        }
    }
    return 0;
}
