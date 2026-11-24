import json, os, pathlib, subprocess, time
import numpy as np

WORKER_ID = 4
SEED_BASE = 202611245000
MAX_SECONDS = 41000
NWORKERS = 4
OUT = pathlib.Path("/kaggle/working")
OUT.mkdir(parents=True, exist_ok=True)

GENERATOR = '#include <stdio.h>\n#include <stdint.h>\n#include <stdlib.h>\n#include <string.h>\n\nstatic int g_bits;\nstatic uint32_t *g_codes;\nstatic int g_count;\n\nstatic void gen_dfs(int pos, int open_n, int int_left, uint32_t cur) {\n    if (pos == g_bits) { g_codes[g_count++] = cur; return; }\n    if (int_left > 0)\n        gen_dfs(pos + 1, open_n + 1, int_left - 1, cur | (1U << (g_bits - 1 - pos)));\n    if (open_n > 1 || (open_n == 1 && int_left == 0))\n        gen_dfs(pos + 1, open_n - 1, int_left, cur);\n}\n\nstatic int cmp_u32(const void *a, const void *b) {\n    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;\n    return (x > y) - (x < y);\n}\n\nstatic inline int bsearch_u32(const uint32_t *arr, int n, uint32_t tgt) {\n    int lo = 0, hi = n - 1;\n    while (lo <= hi) {\n        int mid = lo + ((hi - lo) >> 1);\n        uint32_t v = arr[mid];\n        if (v == tgt) return mid;\n        if (v < tgt) lo = mid + 1; else hi = mid - 1;\n    }\n    return -1;\n}\n\nstatic inline int sub_len(uint32_t code, int bits, int pos) {\n    int bal = 0;\n    for (int p = pos; p < bits; p++) {\n        bal += ((code >> (bits - 1 - p)) & 1U) ? 1 : -1;\n        if (bal < 0) return p - pos + 1;\n    }\n    return bits + 1 - pos;\n}\n\nint main(int argc, char **argv) {\n    int n = (argc > 1) ? atoi(argv[1]) : 16;\n    const char *cnf_path = (argc > 2) ? argv[2] : "A16_sym.cnf";\n    int K = n - 2;\n    g_bits = 2 * K;\n    g_codes = malloc(3000000 * sizeof(uint32_t));\n    g_count = 0;\n    gen_dfs(0, 1, K, 0U);\n    qsort(g_codes, g_count, sizeof(uint32_t), cmp_u32);\n\n    int nv = g_count;\n    int32_t *u = malloc(20000000 * sizeof(int32_t));\n    int32_t *v = malloc(20000000 * sizeof(int32_t));\n    int64_t ne = 0;\n    int bits = g_bits;\n\n    for (int idx = 0; idx < nv; idx++) {\n        uint32_t code = g_codes[idx];\n        for (int p = 0; p < bits; p++) {\n            if (((code >> (bits - 1 - p)) & 1U) == 0U) continue;\n            int lp = p + 1;\n            int ll = sub_len(code, bits, lp);\n            int rp = lp + ll;\n\n            if (lp < bits && ((code >> (bits - 1 - lp)) & 1U) == 1U) {\n                int lenA = sub_len(code, bits, lp + 1);\n                uint32_t nxt = code;\n                for (int i = 0; i < lenA; i++) {\n                    uint32_t b = (code >> (bits - 1 - (lp + 1 + i))) & 1U;\n                    int sh = bits - 1 - (lp + i);\n                    if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);\n                }\n                nxt |= (1U << (bits - 1 - (lp + lenA)));\n                int nb = bsearch_u32(g_codes, nv, nxt);\n                if (idx < nb) { u[ne] = idx; v[ne] = nb; ne++; }\n            }\n            if (rp < bits && ((code >> (bits - 1 - rp)) & 1U) == 1U) {\n                int lenA = ll;\n                uint32_t nxt = code | (1U << (bits - 1 - lp));\n                for (int i = 0; i < lenA; i++) {\n                    uint32_t b = (code >> (bits - 1 - (lp + i))) & 1U;\n                    int sh = bits - 1 - (lp + 1 + i);\n                    if (b) nxt |= (1U << sh); else nxt &= ~(1U << sh);\n                }\n                int nb = bsearch_u32(g_codes, nv, nxt);\n                if (idx < nb) { u[ne] = idx; v[ne] = nb; ne++; }\n            }\n        }\n    }\n    printf("Built A_%d: V=%d, E=%lld\\n", n, nv, (long long)ne);\n\n    FILE *fg = fopen("A16_graph.bin", "wb");\n    fwrite(&nv, sizeof(int32_t), 1, fg);\n    fwrite(&ne, sizeof(int64_t), 1, fg);\n    fwrite(g_codes, sizeof(uint32_t), nv, fg);\n    fwrite(u, sizeof(int32_t), ne, fg);\n    fwrite(v, sizeof(int32_t), ne, fg);\n    fclose(fg);\n\n    /* Break all 24 color symmetries via path v0 - v1 - v2 */\n    int v0 = 0, v1 = v, v2 = -1;\n    for (int64_t e = 0; e < ne; e++) {\n        if (u[e] == v1 && v[e] != v0) { v2 = v[e]; break; }\n        if (v[e] == v1 && u[e] != v0) { v2 = u[e]; break; }\n    }\n\n    FILE *fp = fopen(cnf_path, "w");\n    setvbuf(fp, NULL, _IOFBF, 16 * 1024 * 1024);\n    int64_t nvars = (int64_t)nv * 4;\n    int64_t ncls = (int64_t)nv * 7 + ne * 4 + 3;\n    fprintf(fp, "p cnf %lld %lld\\n", (long long)nvars, (long long)ncls);\n\n    for (int i = 0; i < nv; i++) {\n        int b = 4 * i;\n        fprintf(fp, "%d %d %d %d 0\\n-%d -%d 0\\n-%d -%d 0\\n-%d -%d 0\\n-%d -%d 0\\n-%d -%d 0\\n-%d -%d 0\\n",\n                b+1, b+2, b+3, b+4,\n                b+1, b+2, b+1, b+3, b+1, b+4,\n                b+2, b+3, b+2, b+4, b+3, b+4);\n    }\n    fprintf(fp, "%d 0\\n%d 0\\n-%d 0\\n", 4*v0+1, 4*v1+2, 4*v2+4);\n\n    for (int64_t e = 0; e < ne; e++) {\n        int a = 4 * u[e], b = 4 * v[e];\n        fprintf(fp, "-%d -%d 0\\n-%d -%d 0\\n-%d -%d 0\\n-%d -%d 0\\n",\n                a+1, b+1, a+2, b+2, a+3, b+3, a+4, b+4);\n    }\n    fclose(fp);\n    printf("Wrote %s (%lld vars, %lld clauses)\\n", cnf_path, (long long)nvars, (long long)ncls);\n    return 0;\n}\n'
SOLVER = '/*\n * Direct 4-color stochastic local search for A16_graph.bin.\n * C11, no external libraries. Input format (little-endian host assumed):\n *   int32 nv, int64 ne, uint32 codes[nv], int32 u[ne], int32 v[ne]\n * Output is exactly nv bytes, one color (0..3) per vertex.\n */\n#define _POSIX_C_SOURCE 200809L\n#include <errno.h>\n#include <inttypes.h>\n#include <stdint.h>\n#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>\n#include <time.h>\n\nstatic uint64_t rng_state;\nstatic uint64_t rng(void){ uint64_t x=rng_state; x^=x>>12; x^=x<<25; x^=x>>27; rng_state=x; return x*UINT64_C(2685821657736338717); }\nstatic uint32_t pick(uint32_t n){ return n ? (uint32_t)(rng()%n) : 0; }\nstatic void *xmalloc(size_t n){ void *p=malloc(n); if(!p){fprintf(stderr,"allocation failed (%zu bytes)\\n",n); exit(2);} return p; }\nstatic int read_exact(FILE*f,void*p,size_t z){return fread(p,1,z,f)==z;}\nstatic void usage(const char*p){fprintf(stderr,"Usage: %s GRAPH.bin --seed N [--max-flips N] [--restart-flips N] [--out FILE]\\n",p);}\n\nint main(int argc,char**argv){\n if(argc<2){usage(argv[0]);return 2;}\n const char *graph=argv[1], *out="a16_solution.colors";\n uint64_t seed=1, maxflips=UINT64_MAX, restartflips=10000000; double maxseconds=41000.0; int have_seed=0;\n for(int i=2;i<argc;i++){\n  if(!strcmp(argv[i],"--seed")&&i+1<argc){seed=strtoull(argv[++i],NULL,10);have_seed=1;}\n  else if(!strcmp(argv[i],"--max-flips")&&i+1<argc)maxflips=strtoull(argv[++i],NULL,10);\\n  else if(!strcmp(argv[i],"--max-seconds")&&i+1<argc)maxseconds=strtod(argv[++i],NULL);\n  else if(!strcmp(argv[i],"--restart-flips")&&i+1<argc)restartflips=strtoull(argv[++i],NULL,10);\n  else if(!strcmp(argv[i],"--out")&&i+1<argc)out=argv[++i]; else {usage(argv[0]);return 2;}\n }\n if(!have_seed){fprintf(stderr,"--seed is required\\n");return 2;} rng_state=seed?seed:UINT64_C(0x9e3779b97f4a7c15);\n FILE*f=fopen(graph,"rb"); if(!f){perror(graph);return 2;}\n int32_t nv32; int64_t ne64;\n if(!read_exact(f,&nv32,4)||!read_exact(f,&ne64,8)||nv32<=0||ne64<0){fprintf(stderr,"bad graph header\\n");return 2;}\n uint32_t nv=(uint32_t)nv32; uint64_t ne=(uint64_t)ne64;\n if(ne>SIZE_MAX/(2*sizeof(uint32_t))){fprintf(stderr,"graph too large\\n");return 2;}\n if(fseeko(f,(off_t)((uint64_t)nv*4),SEEK_CUR)){fprintf(stderr,"bad codes section\\n");return 2;}\n uint32_t *eu=xmalloc((size_t)ne*4), *ev=xmalloc((size_t)ne*4);\n for(uint64_t e=0;e<ne;e++){int32_t a;if(!read_exact(f,&a,4)){fprintf(stderr,"truncated u array at %"PRIu64"\\n",e);return 2;} if(a<0||(uint32_t)a>=nv){fprintf(stderr,"invalid u endpoint at %"PRIu64"\\n",e);return 2;}eu[e]=(uint32_t)a;}\n for(uint64_t e=0;e<ne;e++){int32_t b;if(!read_exact(f,&b,4)){fprintf(stderr,"truncated v array at %"PRIu64"\\n",e);return 2;} if(b<0||(uint32_t)b>=nv||eu[e]==(uint32_t)b){fprintf(stderr,"invalid v endpoint at %"PRIu64"\\n",e);return 2;}ev[e]=(uint32_t)b;}\n fclose(f);\n if(ne>UINT32_MAX/2){fprintf(stderr,"too many adjacency entries\\n");return 2;}\n uint32_t *deg=calloc(nv,sizeof(*deg)); if(!deg){perror("calloc");return 2;}\n for(uint64_t e=0;e<ne;e++){deg[eu[e]]++;deg[ev[e]]++;}\n uint64_t *off=xmalloc((size_t)(nv+1)*8); off[0]=0;for(uint32_t v=0;v<nv;v++)off[v+1]=off[v]+deg[v];\n uint32_t *adj=xmalloc((size_t)(2*ne)*4);uint64_t *cur=xmalloc((size_t)nv*8);memcpy(cur,off,(size_t)nv*8);\n for(uint64_t e=0;e<ne;e++){adj[cur[eu[e]]++]=ev[e];adj[cur[ev[e]]++]=eu[e];}free(cur);\n uint8_t *col=xmalloc(nv); uint64_t *tabu=xmalloc((size_t)nv*4*sizeof(uint64_t)); memset(tabu,0,(size_t)nv*4*sizeof(uint64_t));\n uint32_t *bad=xmalloc((size_t)nv*4),*pos=xmalloc((size_t)nv*4),*conf=xmalloc((size_t)nv*4);uint32_t nbad=0;\n for(uint32_t v=0;v<nv;v++){col[v]=(uint8_t)pick(4);pos[v]=UINT32_MAX;}\n uint64_t conflicts=0;\n for(uint32_t v=0;v<nv;v++){uint32_t c=0;for(uint64_t j=off[v];j<off[v+1];j++)c+=(col[adj[j]]==col[v]);conf[v]=c;if(c){bad[nbad]=v;pos[v]=nbad++;}}\n for(uint64_t e=0;e<ne;e++)conflicts+=(col[eu[e]]==col[ev[e]]);\n fprintf(stderr,"graph nv=%u ne=%"PRIu64" seed=%"PRIu64" start_conflicts=%"PRIu64"\\n",nv,ne,seed,conflicts);\n uint64_t steps=0, since_restart=0, best=conflicts; uint32_t stale=0;\n while(conflicts && steps<maxflips){\n  if(nbad==0){fprintf(stderr,"internal bad-list inconsistency\\n");return 3;}\n  uint32_t v=bad[pick(nbad)], old=col[v], counts[4]={0,0,0,0};\n  for(uint64_t j=off[v];j<off[v+1];j++)counts[col[adj[j]]]++;\n  uint32_t bestc=UINT32_MAX, choices[4], nc=0;\n  for(uint32_t c=0;c<4;c++)if(c!=old && (tabu[4*(size_t)v+c] <= steps || counts[c]==0)){\n    if(counts[c]<bestc){bestc=counts[c];nc=0;} if(counts[c]==bestc)choices[nc++]=c;\n  }\n  if(!nc){for(uint32_t c=0;c<4;c++)if(c!=old){if(counts[c]<bestc){bestc=counts[c];nc=0;}if(counts[c]==bestc)choices[nc++]=c;}}\n  uint32_t nw=choices[pick(nc)];\n  /* Tabu is stored modulo 256; short tenure keeps the byte representation exact. */\n  tabu[4*(size_t)v+old]=steps+5+pick(8);\n  col[v]=(uint8_t)nw;\n  for(uint64_t j=off[v];j<off[v+1];j++){\n   uint32_t u=adj[j], before=conf[u];\n   if(col[u]==old) conf[u]--;\n   if(col[u]==nw) conf[u]++;\n   if(!before&&conf[u]){pos[u]=nbad;bad[nbad++]=u;} else if(before&&!conf[u]){uint32_t p=pos[u],last=bad[--nbad];bad[p]=last;pos[last]=p;pos[u]=UINT32_MAX;}\n  }\n  uint32_t cv=0;for(uint64_t j=off[v];j<off[v+1];j++)cv+=(col[adj[j]]==nw);conf[v]=cv;\n  if(cv==0){uint32_t p=pos[v],last=bad[--nbad];bad[p]=last;pos[last]=p;pos[v]=UINT32_MAX;}\n  else { /* v stays conflicted; its list position is unchanged */ }\n  /* Update total conflicts from the incident edges whose endpoint is v. */\n  conflicts = conflicts - counts[old] + counts[nw];\n  steps++; since_restart++;\n  if(steps%1000000==0)fprintf(stderr,"step=%"PRIu64" conflicts=%"PRIu64" best=%"PRIu64"\\n",steps,conflicts,best);\n  if(conflicts<best){best=conflicts;stale=0;}else stale++;\n  if(conflicts==0)break;\n  if(restartflips&&since_restart>=restartflips){for(uint32_t q=0;q<nv;q++){col[q]=(uint8_t)pick(4);tabu[4*(size_t)q]=tabu[4*(size_t)q+1]=tabu[4*(size_t)q+2]=tabu[4*(size_t)q+3]=0;}\n    nbad=0;memset(pos,0xff,(size_t)nv*4);conflicts=0;\n    for(uint32_t q=0;q<nv;q++){uint32_t c=0;for(uint64_t j=off[q];j<off[q+1];j++)c+=(col[adj[j]]==col[q]);conf[q]=c;if(c){pos[q]=nbad;bad[nbad++]=q;}conflicts+=c;}\n    conflicts/=2;since_restart=0;if(conflicts<best)best=conflicts;fprintf(stderr,"restart step=%"PRIu64" conflicts=%"PRIu64"\\n",steps,conflicts);\n  }\n }\n if(conflicts==0){for(uint64_t e=0;e<ne;e++)if(col[eu[e]]==col[ev[e]]){fprintf(stderr,"verification failed at edge %"PRIu64"\\n",e);return 3;}FILE*g=fopen(out,"wb");if(!g){perror(out);return 2;}if(fwrite(col,1,nv,g)!=nv){fprintf(stderr,"write failed\\n");return 2;}fclose(g);printf("SAT seed=%"PRIu64" steps=%"PRIu64" output=%s\\n",seed,steps,out);return 0;}\n fprintf(stderr,"UNKNOWN seed=%"PRIu64" steps=%"PRIu64" conflicts=%"PRIu64" best=%"PRIu64"\\n",seed,steps,conflicts,best);return 1;\n}\n'

def write(path, data):
    pathlib.Path(path).write_text(data)

t0 = time.time()
write("/tmp/build_a16.c", GENERATOR)
write("/tmp/a16_sls.c", SOLVER)
subprocess.run(["gcc", "-O3", "-DNDEBUG", "-std=c11", "/tmp/build_a16.c", "-o", "/tmp/build_a16"], check=True)
subprocess.run(["gcc", "-O3", "-DNDEBUG", "-std=c11", "/tmp/a16_sls.c", "-o", "/tmp/a16_sls"], check=True)
subprocess.run(["/tmp/build_a16", "16", "/tmp/unused.cnf"], check=True)
graph = pathlib.Path("A16_graph.bin")
if not graph.exists():
    raise RuntimeError("graph generator did not create A16_graph.bin")

# One process per available CPU. Individual solver processes share the graph
# but use distinct deterministic seeds and stop at the same wall-clock deadline.
procs = []
logs = []
for i in range(NWORKERS):
    seed = SEED_BASE + i
    raw = OUT / f"A16_colors_seed_{seed}.bin"
    log = open(OUT / f"A16_worker_{i}.log", "w")
    p = subprocess.Popen(
        ["/tmp/a16_sls", str(graph), "--seed", str(seed),
         "--max-seconds", str(MAX_SECONDS), "--out", str(raw)],
        stdout=log, stderr=subprocess.STDOUT
    )
    procs.append((seed, p, raw, log))
    logs.append(log)

winner = None
while procs and time.time() - t0 < MAX_SECONDS:
    for item in list(procs):
        seed, p, raw, log = item
        if p.poll() is not None:
            log.flush()
            if p.returncode == 0 and raw.exists():
                winner = (seed, raw)
                break
            procs.remove(item)
    if winner:
        break
    time.sleep(2)

for seed, p, raw, log in procs:
    if p.poll() is None:
        p.terminate()
for seed, p, raw, log in procs:
    try:
        p.wait(timeout=10)
    except subprocess.TimeoutExpired:
        p.kill()
    log.close()
for log in logs:
    if not log.closed:
        log.close()

status = "UNKNOWN"
verified = False
if winner:
    seed, raw = winner
    colors = np.fromfile(raw, dtype=np.uint8)
    with open(graph, "rb") as f:
        nv = int(np.fromfile(f, dtype=np.int32, count=1)[0])
        ne = int(np.fromfile(f, dtype=np.int64, count=1)[0])
        codes = np.fromfile(f, dtype=np.uint32, count=nv)
        u = np.fromfile(f, dtype=np.int32, count=ne)
        v = np.fromfile(f, dtype=np.int32, count=ne)
    verified = bool(len(colors) == nv and np.all(colors[u] != colors[v]))
    if verified:
        np.savez_compressed(OUT / "A16_4coloring_verified.npz", codes=codes, colors=colors)
        status = "SAT_VERIFIED"
    else:
        status = "ERROR_VERIFICATION_FAILED"

with open(OUT / "A16_kaggle_result.json", "w") as f:
    json.dump({
        "worker_id": WORKER_ID, "seed_base": SEED_BASE,
        "solver": "direct-tabucol-sls", "status": status,
        "verified_4_coloring": verified,
        "elapsed_seconds": round(time.time() - t0, 2),
        "max_seconds": MAX_SECONDS
    }, f, indent=2)
