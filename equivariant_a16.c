#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void code_to_diags(uint32_t code, int n_gon, int diags[][2]) {
    int node_type[64], node_l[64], node_r[64];
    int nsp = 0, leaf = 0, d_idx = 0;
    int bits = 2 * (n_gon - 2);

    for (int p = 0; p < bits; p++) {
        int bit = (int)((code >> (bits - 1 - p)) & 1U);
        node_type[nsp] = bit;
        if (bit == 0) {
            node_l[nsp] = leaf;
            node_r[nsp] = leaf + 1;
            leaf++;
        }
        nsp++;
        while (nsp >= 3 && node_type[nsp - 3] == 1 && node_type[nsp - 2] == 0 && node_type[nsp - 1] == 0) {
            int l = node_l[nsp - 2], r = node_r[nsp - 1];
            if (!(l == 0 && r == n_gon - 1)) {
                diags[d_idx][0] = l;
                diags[d_idx][1] = r;
                d_idx++;
            }
            nsp -= 2;
            node_type[nsp - 1] = 0;
            node_l[nsp - 1] = l;
            node_r[nsp - 1] = r;
        }
    }
    node_type[nsp] = 0;
    node_l[nsp] = leaf;
    node_r[nsp] = leaf + 1;
    nsp++;
    while (nsp >= 3 && node_type[nsp - 3] == 1 && node_type[nsp - 2] == 0 && node_type[nsp - 1] == 0) {
        int l = node_l[nsp - 2], r = node_r[nsp - 1];
        if (!(l == 0 && r == n_gon - 1)) {
            diags[d_idx][0] = l;
            diags[d_idx][1] = r;
            d_idx++;
        }
        nsp -= 2;
        node_type[nsp - 1] = 0;
        node_l[nsp - 1] = l;
        node_r[nsp - 1] = r;
    }
}

static uint32_t diags_to_code(int n_gon, int diags[][2]) {
    uint8_t adj[32][32];
    memset(adj, 0, sizeof(adj));
    for (int i = 0; i < n_gon - 3; i++) {
        int a = diags[i][0], b = diags[i][1];
        adj[a][b] = adj[b][a] = 1;
    }
    uint32_t code = 0;
    int pos = 0, bits = 2 * (n_gon - 2);
    int sl[64], sr[64], sp = 0;
    sl[sp] = 0; sr[sp] = n_gon - 1; sp++;

    while (sp > 0) {
        sp--;
        int l = sl[sp], r = sr[sp];
        if (r - l == 1) {
            if (pos < bits) pos++;
        } else {
            int mid = -1;
            for (int k = l + 1; k < r; k++) {
                int ok_l = (k - l == 1) || adj[l][k];
                int ok_r = (r - k == 1) || adj[k][r];
                if (ok_l && ok_r) { mid = k; break; }
            }
            if (mid < 0) {
                fprintf(stderr, "ERROR: invalid triangulation interval [%d, %d]\n", l, r);
                exit(1);
            }
            code |= (1U << (bits - 1 - pos));
            pos++;
            sl[sp] = mid; sr[sp] = r; sp++;
            sl[sp] = l;   sr[sp] = mid; sp++;
        }
    }
    return code;
}

static int bsearch_u32(const uint32_t *arr, int n, uint32_t tgt) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (arr[mid] == tgt) return mid;
        if (arr[mid] < tgt) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <verify_only|build|lift> ...\n", argv[0]);
        return 1;
    }

    FILE *fg = fopen("A16_graph.bin", "rb");
    if (!fg) { perror("A16_graph.bin"); return 1; }
    int32_t nv; int64_t ne;
    if (fread(&nv, sizeof(int32_t), 1, fg) != 1 || fread(&ne, sizeof(int64_t), 1, fg) != 1) {
        fprintf(stderr, "Failed reading header\n"); return 1;
    }
    uint32_t *codes = malloc((size_t)nv * sizeof(uint32_t));
    int32_t *eu = malloc((size_t)ne * sizeof(int32_t));
    int32_t *ev = malloc((size_t)ne * sizeof(int32_t));
    if (fread(codes, sizeof(uint32_t), (size_t)nv, fg) != (size_t)nv ||
        fread(eu, sizeof(int32_t), (size_t)ne, fg) != (size_t)ne ||
        fread(ev, sizeof(int32_t), (size_t)ne, fg) != (size_t)ne) {
        fprintf(stderr, "Failed reading graph arrays\n"); return 1;
    }
    fclose(fg);

    int32_t *rot1 = malloc((size_t)nv * sizeof(int32_t));
    uint8_t *seen = calloc((size_t)nv, 1);
    for (int i = 0; i < nv; i++) {
        int diags[16][2];
        code_to_diags(codes[i], 16, diags);
        if (diags_to_code(16, diags) != codes[i]) {
            fprintf(stderr, "MATH ERROR: round-trip code mismatch at vertex %d\n", i);
            return 1;
        }
        for (int d = 0; d < 13; d++) {
            int a = (diags[d][0] + 1) & 15;
            int b = (diags[d][1] + 1) & 15;
            if (a > b) { int t = a; a = b; b = t; }
            diags[d][0] = a; diags[d][1] = b;
        }
        uint32_t nc = diags_to_code(16, diags);
        int idx = bsearch_u32(codes, nv, nc);
        if (idx < 0 || seen[idx]) {
            fprintf(stderr, "MATH ERROR: rot1 is not a bijection at vertex %d\n", i);
            return 1;
        }
        seen[idx] = 1;
        rot1[i] = idx;
    }
    free(seen);

    for (int i = 0; i < nv; i++) {
        int cur = i;
        for (int s = 0; s < 16; s++) cur = rot1[cur];
        if (cur != i) {
            fprintf(stderr, "MATH ERROR: rot1^16 != id at vertex %d\n", i);
            return 1;
        }
    }

    uint64_t *edge_set = malloc((size_t)ne * sizeof(uint64_t));
    for (int64_t e = 0; e < ne; e++) {
        uint64_t a = (uint64_t)(eu[e] < ev[e] ? eu[e] : ev[e]);
        uint64_t b = (uint64_t)(eu[e] < ev[e] ? ev[e] : eu[e]);
        edge_set[e] = (a << 32) | b;
    }
    qsort(edge_set, (size_t)ne, sizeof(uint64_t), cmp_u64);
    for (int64_t e = 0; e < ne; e++) {
        int ru = rot1[eu[e]], rv = rot1[ev[e]];
        uint64_t key = ((uint64_t)(ru < rv ? ru : rv) << 32) | (uint64_t)(ru < rv ? rv : ru);
        int lo = 0, hi = (int)ne - 1, found = 0;
        while (lo <= hi) {
            int mid = lo + ((hi - lo) >> 1);
            if (edge_set[mid] == key) { found = 1; break; }
            if (edge_set[mid] < key) lo = mid + 1; else hi = mid - 1;
        }
        if (!found) {
            fprintf(stderr, "MATH ERROR: edge (%d,%d) mapped to non-edge (%d,%d)\n", eu[e], ev[e], ru, rv);
            return 1;
        }
    }
    free(edge_set);

    if (strcmp(argv[1], "verify_only") == 0) {
        printf("VERIFIED: rot1 is an exact order-16 graph automorphism of A16 (nv=%d, ne=%lld).\n", nv, (long long)ne);
        return 0;
    }

    if (argc < 4) {
        fprintf(stderr, "Usage: %s <build|lift> <k_rot> <shift> <file1> [file2]\n", argv[0]);
        return 1;
    }

    int k_rot = atoi(argv[2]);
    int shift = atoi(argv[3]) & 3;
    int step = 16 / k_rot;

    int32_t *rho = malloc((size_t)nv * sizeof(int32_t));
    for (int i = 0; i < nv; i++) {
        int cur = i;
        for (int s = 0; s < step; s++) cur = rot1[cur];
        rho[i] = cur;
    }

    int32_t *orb_id = malloc((size_t)nv * sizeof(int32_t));
    int8_t  *orb_pos = malloc((size_t)nv * sizeof(int8_t));
    memset(orb_id, -1, (size_t)nv * sizeof(int32_t));
    int num_orbs = 0;
    for (int i = 0; i < nv; i++) {
        if (orb_id[i] != -1) continue;
        int cur = i, len = 0;
        while (orb_id[cur] == -1) {
            orb_id[cur] = num_orbs;
            orb_pos[cur] = (int8_t)len;
            len++;
            cur = rho[cur];
        }
        if (((len * shift) & 3) != 0) {
            printf("REJECTED (k=%2d, s=%d): orbit %d has len=%d, len*s=%d != 0 mod 4\n", k_rot, shift, num_orbs, len, len * shift);
            return 2;
        }
        num_orbs++;
    }

    if (strcmp(argv[1], "lift") == 0) {
        if (argc < 6) {
            fprintf(stderr, "Usage: %s lift <k_rot> <shift> <quotient_log> <out_full_log>\n", argv[0]);
            return 1;
        }
        const char *qlog = argv[4];
        const char *out_full = argv[5];
        int8_t *orb_col = malloc((size_t)num_orbs);
        memset(orb_col, -1, (size_t)num_orbs);

        FILE *fl = fopen(qlog, "r");
        if (!fl) { perror(qlog); return 1; }
        char tok[64];
        int sat = 0;
        while (fscanf(fl, "%63s", tok) == 1) {
            if (strcmp(tok, "SATISFIABLE") == 0) sat = 1;
            else if (tok[0] == 'v' && tok[1] == '\0') continue;
            else if (tok[0] == 's' && tok[1] == '\0') continue;
            else {
                long lit = strtol(tok, NULL, 10);
                if (lit > 0) {
                    int idx = (int)(lit - 1);
                    int o = idx / 4, c = idx & 3;
                    if (o >= 0 && o < num_orbs && orb_col[o] == -1) orb_col[o] = (int8_t)c;
                }
            }
        }
        fclose(fl);
        if (!sat) { fprintf(stderr, "Log %s is not SATISFIABLE\n", qlog); return 1; }

        int8_t *vcol = malloc((size_t)nv);
        for (int i = 0; i < nv; i++) {
            vcol[i] = (int8_t)((orb_col[orb_id[i]] + orb_pos[i] * shift) & 3);
        }
        for (int64_t e = 0; e < ne; e++) {
            if (vcol[eu[e]] == vcol[ev[e]]) {
                fprintf(stderr, "LIFT VERIFY FAILED on edge %lld: (%d, %d) both color %d\n",
                        (long long)e, eu[e], ev[e], vcol[eu[e]]);
                return 1;
            }
        }
        FILE *fo = fopen(out_full, "w");
        fprintf(fo, "s SATISFIABLE\nv");
        for (int i = 0; i < nv; i++) {
            int b = 4 * i;
            for (int c = 0; c < 4; c++) {
                fprintf(fo, " %d", (vcol[i] == c) ? (b + c + 1) : -(b + c + 1));
            }
            if ((i & 1023) == 1023) fprintf(fo, "\nv");
        }
        fprintf(fo, " 0\n");
        fclose(fo);
        printf("LIFT VERIFIED: All %lld edges valid! Wrote full solution to %s\n", (long long)ne, out_full);
        return 0;
    }

    uint64_t *qedges = malloc((size_t)ne * sizeof(uint64_t));
    int64_t qcount = 0;
    for (int64_t e = 0; e < ne; e++) {
        int u = eu[e], v = ev[e];
        int oa = orb_id[u], ob = orb_id[v];
        int pa = orb_pos[u], pb = orb_pos[v];
        int delta = (((pa - pb) * shift) % 4 + 4) % 4;
        if (oa == ob) {
            if (delta == 0) {
                printf("REJECTED (k=%2d, s=%d): self-loop conflict in orbit %d (edge %d-%d)\n", k_rot, shift, oa, u, v);
                return 3;
            }
            continue;
        }
        if (oa > ob) {
            int t = oa; oa = ob; ob = t;
            delta = (4 - delta) & 3;
        }
        qedges[qcount++] = ((uint64_t)oa << 24) | ((uint64_t)ob << 2) | (uint64_t)delta;
    }
    qsort(qedges, (size_t)qcount, sizeof(uint64_t), cmp_u64);
    int64_t uq = 0;
    for (int64_t i = 0; i < qcount; i++) {
        if (i == 0 || qedges[i] != qedges[i - 1]) qedges[uq++] = qedges[i];
    }

    const char *out_cnf = argv[4];
    FILE *fc = fopen(out_cnf, "w");
    setvbuf(fc, NULL, _IOFBF, 16 * 1024 * 1024);
    fprintf(fc, "p cnf %lld %lld\n", (long long)num_orbs * 4, (long long)num_orbs + uq * 4 + 1);
    for (int i = 0; i < num_orbs; i++) {
        int b = 4 * i;
        fprintf(fc, "%d %d %d %d 0\n", b + 1, b + 2, b + 3, b + 4);
    }
    fprintf(fc, "1 0\n");
    for (int64_t i = 0; i < uq; i++) {
        int oa = (int)(qedges[i] >> 24);
        int ob = (int)((qedges[i] >> 2) & 0x3FFFFF);
        int delta = (int)(qedges[i] & 3);
        for (int c = 0; c < 4; c++) {
            int cb = (c + delta) & 3;
            fprintf(fc, "-%d -%d 0\n", 4 * oa + c + 1, 4 * ob + cb + 1);
        }
    }
    fclose(fc);
    printf("VALID    (k=%2d, s=%d): %d orbits, %lld quotient constraints -> %s\n",
           k_rot, shift, num_orbs, (long long)uq, out_cnf);
    return 0;
}
