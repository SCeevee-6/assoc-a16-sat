import os, json, subprocess, time
import numpy as np

t0 = time.time()
K_ROT = 8
SHIFT = 1
SEEDS = [910001, 910002, 910003, 910004]

subprocess.run(["apt-get", "update", "-qq"], check=True)
subprocess.run(["apt-get", "install", "-y", "-qq", "gcc", "g++", "make", "git"], check=True)

if not os.path.exists("/tmp/kissat/build/kissat"):
    subprocess.run(["git", "clone", "--depth=1", "https://github.com/arminbiere/kissat.git", "/tmp/kissat"], check=True)
    subprocess.run(["./configure", "--quiet"], cwd="/tmp/kissat", check=True)
    subprocess.run(["make", "-j4"], cwd="/tmp/kissat", check=True)

with open("/tmp/build_quot.c", "w") as f:
    f.write(r"""
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int g_bits, g_count;
static uint32_t *g_codes;

static void gen_dfs(int pos, int open_n, int int_left, uint32_t cur) {
    if (pos == g_bits) { g_codes[g_count++] = cur; return; }
    if (int_left > 0)
        gen_dfs(pos + 1, open_n + 1, int_left - 1, cur | (1U << (g_bits - 1 - pos)));
    if (open_n > 1 || (open_n == 1 && int_left == 0))
        gen_dfs(pos + 1, open_n - 1, int_left, cur);
}
static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}
static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}
static inline int bsearch_u32(const uint32_t *arr, int n, uint32_t tgt) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (arr[mid] == tgt) return mid;
        if (arr[mid] < tgt) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}
static inline int sub_len(uint32_t code, int bits, int pos) {
    int bal = 0;
    for (int p = pos; p < bits; p++) {
        bal += ((code >> (bits - 1 - p)) & 1U) ? 1 : -1;
        if (bal < 0) return p - pos + 1;
    }
    return bits + 1 - pos;
}
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
    uint8_t adj[20][20];
    memset(adj, 0, sizeof(adj));
    for (int i = 0; i < n_gon - 3; i++) {
        int a = diags[i][0], b = diags[i][1];
        adj[a][b] = adj[b][a] = 1;
    }
    uint32_t code = 0;
    int pos = 0, bits = 2 * (n_gon - 2), sl[64], sr[64], sp = 0;
    sl[sp] = 0; sr[sp] = n_gon - 1; sp++;
    while (sp > 0) {
        sp--;
        int l = sl[sp], r = sr[sp];
        if (r - l == 1) { if (pos < bits) pos++; }
        else {
            int mid = -1;
            for (int k = l + 1; k < r; k++) {
                if (((k - l == 1) || adj[l][k]) && ((r - k == 1) || adj[k][r])) { mid = k; break; }
            }
            code |= (1U << (bits - 1 - pos)); pos++;
            sl[sp] = mid; sr[sp] = r; sp++;
            sl[sp] = l;   sr[sp] = mid; sp++;
        }
    }
    return code;
}

int main(int argc, char **argv) {
    int k_rot = atoi(argv[1]), shift = atoi(argv[2]) & 3, step = 16 / k_rot;
    g_bits = 28;
    g_codes = malloc(3000000 * sizeof(uint32_t));
    g_count = 0;
    gen_dfs(0, 1, 14, 0U);
    qsort(g_codes, g_count, sizeof(uint32_t), cmp_u32);
    int nv = g_count;

    int32_t *eu = malloc(20000000 * sizeof(int32_t));
    int32_t *ev = malloc(20000000 * sizeof(int32_t));
    int64_t ne = 0;
    for (int idx = 0; idx < nv; idx++) {
        uint32_t code = g_codes[idx];
        for (int p = 0; p < 28; p++) {
            if (((code >> (27 - p)) & 1U) == 0U) continue;
            int lp = p + 1, ll = sub_len(code, 28, lp), rp = lp + ll;
            if (lp < 28 && ((code >> (27 - lp)) & 1U) == 1U) {
                int lenA = sub_len(code, 28, lp + 1);
                uint32_t nxt = code;
                for (int i = 0; i < lenA; i++) {
                    uint32_t b = (code >> (27 - (lp + 1 + i))) & 1U;
                    int sh = 27 - (lp + i);
                    if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);
                }
                nxt |= (1U << (27 - (lp + lenA)));
                int nb = bsearch_u32(g_codes, nv, nxt);
                if (idx < nb) { eu[ne] = idx; ev[ne] = nb; ne++; }
            }
            if (rp < 28 && ((code >> (27 - rp)) & 1U) == 1U) {
                int lenA = ll;
                uint32_t nxt = code | (1U << (27 - lp));
                for (int i = 0; i < lenA; i++) {
                    uint32_t b = (code >> (27 - (lp + i))) & 1U;
                    int sh = 27 - (lp + 1 + i);
                    if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);
                }
                int nb = bsearch_u32(g_codes, nv, nxt);
                if (idx < nb) { eu[ne] = idx; ev[ne] = nb; ne++; }
            }
        }
    }

    int32_t *rho = malloc(nv * sizeof(int32_t));
    for (int i = 0; i < nv; i++) {
        int diags[16][2];
        code_to_diags(g_codes[i], 16, diags);
        for (int d = 0; d < 13; d++) {
            int a = (diags[d][0] + step) & 15, b = (diags[d][1] + step) & 15;
            if (a > b) { int t = a; a = b; b = t; }
            diags[d][0] = a; diags[d][1] = b;
        }
        rho[i] = bsearch_u32(g_codes, nv, diags_to_code(16, diags));
    }

    int32_t *orb_id = malloc(nv * sizeof(int32_t));
    int8_t  *orb_pos = malloc(nv * sizeof(int8_t));
    memset(orb_id, -1, nv * sizeof(int32_t));
    int num_orbs = 0;
    for (int i = 0; i < nv; i++) {
        if (orb_id[i] != -1) continue;
        int cur = i, len = 0;
        while (orb_id[cur] == -1) {
            orb_id[cur] = num_orbs; orb_pos[cur] = (int8_t)len; len++; cur = rho[cur];
        }
        num_orbs++;
    }

    if (argc >= 5 && strcmp(argv[3], "lift") == 0) {
        int8_t *orb_col = malloc(num_orbs);
        memset(orb_col, -1, num_orbs);
        FILE *fl = fopen(argv[4], "r");
        char tok[64]; int sat = 0;
        while (fscanf(fl, "%63s", tok) == 1) {
            if (strcmp(tok, "SATISFIABLE") == 0) sat = 1;
            else {
                long lit = strtol(tok, NULL, 10);
                if (lit > 0) {
                    int idx = (int)(lit - 1), o = idx / 4, c = idx & 3;
                    if (o >= 0 && o < num_orbs && orb_col[o] == -1) orb_col[o] = (int8_t)c;
                }
            }
        }
        fclose(fl);
        if (!sat) return 2;
        int8_t *vcol = malloc(nv);
        for (int i = 0; i < nv; i++) vcol[i] = (int8_t)((orb_col[orb_id[i]] + orb_pos[i] * shift) & 3);
        for (int64_t e = 0; e < ne; e++) {
            if (vcol[eu[e]] == vcol[ev[e]]) return 3;
        }
        FILE *fo = fopen("/kaggle/working/A16_4coloring_verified.bin", "wb");
        fwrite(vcol, 1, nv, fo);
        fclose(fo);
        return 0;
    }

    uint64_t *qedges = malloc(ne * sizeof(uint64_t));
    int64_t qcount = 0;
    for (int64_t e = 0; e < ne; e++) {
        int u = eu[e], v = ev[e], oa = orb_id[u], ob = orb_id[v];
        int delta = (((orb_pos[u] - orb_pos[v]) * shift) % 4 + 4) % 4;
        if (oa == ob) continue;
        if (oa > ob) { int t = oa; oa = ob; ob = t; delta = (4 - delta) & 3; }
        qedges[qcount++] = ((uint64_t)oa << 24) | ((uint64_t)ob << 2) | (uint64_t)delta;
    }
    qsort(qedges, qcount, sizeof(uint64_t), cmp_u64);
    int64_t uq = 0;
    for (int64_t i = 0; i < qcount; i++) {
        if (i == 0 || qedges[i] != qedges[i - 1]) qedges[uq++] = qedges[i];
    }

    FILE *fc = fopen("/tmp/quot.cnf", "w");
    setvbuf(fc, NULL, _IOFBF, 16 * 1024 * 1024);
    fprintf(fc, "p cnf %lld %lld\n", (long long)num_orbs * 4, (long long)num_orbs + uq * 4 + 1);
    for (int i = 0; i < num_orbs; i++) {
        int b = 4 * i;
        fprintf(fc, "%d %d %d %d 0\n", b + 1, b + 2, b + 3, b + 4);
    }
    fprintf(fc, "1 0\n");
    for (int64_t i = 0; i < uq; i++) {
        int oa = (int)(qedges[i] >> 24), ob = (int)((qedges[i] >> 2) & 0x3FFFFF), delta = (int)(qedges[i] & 3);
        for (int c = 0; c < 4; c++) {
            fprintf(fc, "-%d -%d 0\n", 4 * oa + c + 1, 4 * ob + ((c + delta) & 3) + 1);
        }
    }
    fclose(fc);
    return 0;
}
""")

subprocess.run(["gcc", "-O3", "/tmp/build_quot.c", "-o", "/tmp/build_quot"], check=True)
subprocess.run(["/tmp/build_quot", str(K_ROT), str(SHIFT)], check=True)

procs = []
for s in SEEDS:
    logp = f"/tmp/kissat_{s}.log"
    lf = open(logp, "w")
    p = subprocess.Popen(
        ["/tmp/kissat/build/kissat", f"--seed={s}", "--walkinitially", "--time=41000", "/tmp/quot.cnf"],
        stdout=lf, stderr=subprocess.STDOUT
    )
    procs.append((s, p, lf, logp))

winner_seed = None
final_status = "UNKNOWN"
while procs:
    for item in list(procs):
        s, p, lf, logp = item
        if p.poll() is not None:
            lf.close()
            procs.remove(item)
            with open(logp) as rf:
                for line in rf:
                    if line.startswith("s "):
                        st = line.strip()
                        if "SATISFIABLE" in st and "UNSATISFIABLE" not in st:
                            winner_seed = s
                            final_status = st
                            for _, op, olf, _ in procs:
                                op.terminate()
                                olf.close()
                            procs.clear()
                        elif final_status == "UNKNOWN":
                            final_status = st
                        break
    if procs:
        time.sleep(5)

lift_verified = False
if winner_seed is not None:
    rc = subprocess.run(["/tmp/build_quot", str(K_ROT), str(SHIFT), "lift", f"/tmp/kissat_{winner_seed}.log"]).returncode
    if rc == 0:
        lift_verified = True
        arr = np.fromfile("/kaggle/working/A16_4coloring_verified.bin", dtype=np.int8)
        np.savez_compressed("/kaggle/working/A16_4coloring_verified.npz", colors=arr)

with open("/kaggle/working/A16_kaggle_result.json", "w") as jf:
    json.dump({
        "k_rot": K_ROT,
        "shift": SHIFT,
        "seeds": SEEDS,
        "winner_seed": winner_seed,
        "status": final_status,
        "lift_verified": lift_verified,
        "elapsed_sec": round(time.time() - t0, 2)
    }, jf, indent=2)
