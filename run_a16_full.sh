#!/usr/bin/env bash
set -euo pipefail

echo "=== [1/5] Installing build dependencies ==="
sudo apt-get update -qq
sudo apt-get install -y -qq build-essential git python3 python3-numpy

echo "=== [2/5] Writing & Running C11 Associahedron A_16 Generator ==="
cat << 'CEOF' > build_a16.c
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int g_bits;
static uint32_t *g_codes;
static int g_count;

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

static inline int bsearch_u32(const uint32_t *arr, int n, uint32_t tgt) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = lo + ((hi - lo) >> 1);
        uint32_t v = arr[mid];
        if (v == tgt) return mid;
        if (v < tgt) lo = mid + 1; else hi = mid - 1;
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

int main(int argc, char **argv) {
    int n = (argc > 1) ? atoi(argv) : 16;
    const char *cnf_path = (argc > 2) ? argv : "A16_sym.cnf";
    int K = n - 2;
    g_bits = 2 * K;
    g_codes = malloc(3000000 * sizeof(uint32_t));
    g_count = 0;
    gen_dfs(0, 1, K, 0U);
    qsort(g_codes, g_count, sizeof(uint32_t), cmp_u32);

    int nv = g_count;
    int32_t *u = malloc(20000000 * sizeof(int32_t));
    int32_t *v = malloc(20000000 * sizeof(int32_t));
    int64_t ne = 0;
    int bits = g_bits;

    for (int idx = 0; idx < nv; idx++) {
        uint32_t code = g_codes[idx];
        for (int p = 0; p < bits; p++) {
            if (((code >> (bits - 1 - p)) & 1U) == 0U) continue;
            int lp = p + 1;
            int ll = sub_len(code, bits, lp);
            int rp = lp + ll;

            if (lp < bits && ((code >> (bits - 1 - lp)) & 1U) == 1U) {
                int lenA = sub_len(code, bits, lp + 1);
                uint32_t nxt = code;
                for (int i = 0; i < lenA; i++) {
                    uint32_t b = (code >> (bits - 1 - (lp + 1 + i))) & 1U;
                    int sh = bits - 1 - (lp + i);
                    if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);
                }
                nxt |= (1U << (bits - 1 - (lp + lenA)));
                int nb = bsearch_u32(g_codes, nv, nxt);
                if (idx < nb) { u[ne] = idx; v[ne] = nb; ne++; }
            }
            if (rp < bits && ((code >> (bits - 1 - rp)) & 1U) == 1U) {
                int lenA = ll;
                uint32_t nxt = code | (1U << (bits - 1 - lp));
                for (int i = 0; i < lenA; i++) {
                    uint32_t b = (code >> (bits - 1 - (lp + i))) & 1U;
                    int sh = bits - 1 - (lp + 1 + i);
                    if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);
                }
                int nb = bsearch_u32(g_codes, nv, nxt);
                if (idx < nb) { u[ne] = idx; v[ne] = nb; ne++; }
            }
        }
    }
    printf("Built A_%d: V=%d, E=%lld\n", n, nv, (long long)ne);

    FILE *fg = fopen("A16_graph.bin", "wb");
    fwrite(&nv, sizeof(int32_t), 1, fg);
    fwrite(&ne, sizeof(int64_t), 1, fg);
    fwrite(g_codes, sizeof(uint32_t), nv, fg);
    fwrite(u, sizeof(int32_t), ne, fg);
    fwrite(v, sizeof(int32_t), ne, fg);
    fclose(fg);

    /* Break all 24 color symmetries via path v0 - v1 - v2 */
    int v0 = 0, v1 = v, v2 = -1;
    for (int64_t e = 0; e < ne; e++) {
        if (u[e] == v1 && v[e] != v0) { v2 = v[e]; break; }
        if (v[e] == v1 && u[e] != v0) { v2 = u[e]; break; }
    }

    FILE *fp = fopen(cnf_path, "w");
    setvbuf(fp, NULL, _IOFBF, 16 * 1024 * 1024);
    int64_t nvars = (int64_t)nv * 4;
    int64_t ncls = (int64_t)nv * 7 + ne * 4 + 3;
    fprintf(fp, "p cnf %lld %lld\n", (long long)nvars, (long long)ncls);

    for (int i = 0; i < nv; i++) {
        int b = 4 * i;
        fprintf(fp, "%d %d %d %d 0\n-%d -%d 0\n-%d -%d 0\n-%d -%d 0\n-%d -%d 0\n-%d -%d 0\n-%d -%d 0\n",
                b+1, b+2, b+3, b+4,
                b+1, b+2, b+1, b+3, b+1, b+4,
                b+2, b+3, b+2, b+4, b+3, b+4);
    }
    fprintf(fp, "%d 0\n%d 0\n-%d 0\n", 4*v0+1, 4*v1+2, 4*v2+4);

    for (int64_t e = 0; e < ne; e++) {
        int a = 4 * u[e], b = 4 * v[e];
        fprintf(fp, "-%d -%d 0\n-%d -%d 0\n-%d -%d 0\n-%d -%d 0\n",
                a+1, b+1, a+2, b+2, a+3, b+3, a+4, b+4);
    }
    fclose(fp);
    printf("Wrote %s (%lld vars, %lld clauses)\n", cnf_path, (long long)nvars, (long long)ncls);
    return 0;
}
CEOF

if [ ! -f "A16_sym.cnf" ]; then
    gcc -O3 -march=native build_a16.c -o build_a16
    ./build_a16 16 A16_sym.cnf
else
    echo "A16_sym.cnf already exists on disk, skipping rebuild."
fi

echo "=== [3/5] Cloning & Compiling Plingeling (4-core shared-memory SAT solver) ==="
GH_USER="arminbiere"
GH_REPO="lingeling"
if [ ! -d "lingeling" ]; then
    git clone -q "https://github.com{GH_USER}/${GH_REPO}.git"
    cd lingeling
    ./configure.sh
    make -j$(nproc) plingeling
    cd ..
fi

echo "=== [4/5] Running Plingeling across $(nproc) cores ==="
START_TS=$(date +%s)
set +e
./lingeling/plingeling -t $(nproc) A16_sym.cnf | tee plingeling_a16.log
EXIT_CODE=$?
set -e
END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

echo "=== [5/5] Parsing & Verifying Output (Exit code: $EXIT_CODE, Elapsed: ${ELAPSED}s) ==="
python3 - << PYEOF
import json, struct
import numpy as np

with open("A16_graph.bin", "rb") as f:
    nv = struct.unpack("<i", f.read(4))
    ne = struct.unpack("<q", f.read(8))
    codes = np.frombuffer(f.read(nv * 4), dtype=np.uint32)
    u = np.frombuffer(f.read(ne * 4), dtype=np.int32)
    v = np.frombuffer(f.read(ne * 4), dtype=np.int32)

status = "UNKNOWN"
pos_lits = []
with open("plingeling_a16.log", "r") as f:
    for line in f:
        if line.startswith("s "):
            status = line.strip()
        elif line.startswith("v "):
            for tok in line.split()[1:]:
                val = int(tok)
                if val > 0:
                    pos_lits.append(val - 1)

verified = False
if "SATISFIABLE" in status and "UNSATISFIABLE" not in status:
    pos = np.array(pos_lits, dtype=np.int32)
    colors = np.full(nv, -1, dtype=np.int8)
    colors[pos // 4] = (pos % 4).astype(np.int8)
    verified = bool(np.all((colors >= 0) & (colors <= 3)) and np.all(colors[u] != colors[v]))
    if verified:
        np.savez_compressed("A16_4coloring_verified.npz", codes=codes, colors=colors)
        print("SUCCESS: CHI(A_16) = 4 VERIFIED! Saved A16_4coloring_verified.npz")
    else:
        print("ERROR: Model returned SATISFIABLE but failed verification!")
elif "UNSATISFIABLE" in status:
    print("SUCCESS: CHI(A_16) >= 5 PROVEN (UNSATISFIABLE)!")
else:
    print("Solver finished with status:", status)

with open("A16_result.json", "w") as f:
    json.dump({
        "n": 16,
        "V": int(nv),
        "E": int(ne),
        "solver": "plingeling",
        "threads": $(nproc),
        "status": status,
        "verified_4_coloring": verified,
        "solve_seconds": $ELAPSED
    }, f, indent=2)
PYEOF
