/*
 * Direct 4-color stochastic local search for A16_graph.bin.
 * C11, no external libraries. Input format (little-endian host assumed):
 *   int32 nv, int64 ne, uint32 codes[nv], int32 u[ne], int32 v[ne]
 * Output is exactly nv bytes, one color (0..3) per vertex.
 */
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static uint64_t rng_state;
static uint64_t rng(void){ uint64_t x=rng_state; x^=x>>12; x^=x<<25; x^=x>>27; rng_state=x; return x*UINT64_C(2685821657736338717); }
static uint32_t pick(uint32_t n){ return n ? (uint32_t)(rng()%n) : 0; }
static void *xmalloc(size_t n){ void *p=malloc(n); if(!p){fprintf(stderr,"allocation failed (%zu bytes)\n",n); exit(2);} return p; }
static int read_exact(FILE*f,void*p,size_t z){return fread(p,1,z,f)==z;}
static void usage(const char*p){fprintf(stderr,"Usage: %s GRAPH.bin --seed N [--max-flips N] [--restart-flips N] [--out FILE]\n",p);}

int main(int argc,char**argv){
 if(argc<2){usage(argv[0]);return 2;}
 const char *graph=argv[1], *out="a16_solution.colors";
 uint64_t seed=1, maxflips=UINT64_MAX, restartflips=10000000; double maxseconds=41000.0; int have_seed=0;
 for(int i=2;i<argc;i++){
  if(!strcmp(argv[i],"--seed")&&i+1<argc){seed=strtoull(argv[++i],NULL,10);have_seed=1;}
  else if(!strcmp(argv[i],"--max-flips")&&i+1<argc)maxflips=strtoull(argv[++i],NULL,10);
  else if(!strcmp(argv[i],"--max-seconds")&&i+1<argc)maxseconds=strtod(argv[++i],NULL);
  else if(!strcmp(argv[i],"--restart-flips")&&i+1<argc)restartflips=strtoull(argv[++i],NULL,10);
  else if(!strcmp(argv[i],"--out")&&i+1<argc)out=argv[++i]; else {usage(argv[0]);return 2;}
 }
 if(!have_seed){fprintf(stderr,"--seed is required\n");return 2;} rng_state=seed?seed:UINT64_C(0x9e3779b97f4a7c15);
 FILE*f=fopen(graph,"rb"); if(!f){perror(graph);return 2;}
 int32_t nv32; int64_t ne64;
 if(!read_exact(f,&nv32,4)||!read_exact(f,&ne64,8)||nv32<=0||ne64<0){fprintf(stderr,"bad graph header\n");return 2;}
 uint32_t nv=(uint32_t)nv32; uint64_t ne=(uint64_t)ne64;
 if(ne>SIZE_MAX/(2*sizeof(uint32_t))){fprintf(stderr,"graph too large\n");return 2;}
 if(fseeko(f,(off_t)((uint64_t)nv*4),SEEK_CUR)){fprintf(stderr,"bad codes section\n");return 2;}
 uint32_t *eu=xmalloc((size_t)ne*4), *ev=xmalloc((size_t)ne*4);
 for(uint64_t e=0;e<ne;e++){int32_t a;if(!read_exact(f,&a,4)){fprintf(stderr,"truncated u array at %"PRIu64"\n",e);return 2;} if(a<0||(uint32_t)a>=nv){fprintf(stderr,"invalid u endpoint at %"PRIu64"\n",e);return 2;}eu[e]=(uint32_t)a;}
 for(uint64_t e=0;e<ne;e++){int32_t b;if(!read_exact(f,&b,4)){fprintf(stderr,"truncated v array at %"PRIu64"\n",e);return 2;} if(b<0||(uint32_t)b>=nv||eu[e]==(uint32_t)b){fprintf(stderr,"invalid v endpoint at %"PRIu64"\n",e);return 2;}ev[e]=(uint32_t)b;}
 fclose(f);
 if(ne>UINT32_MAX/2){fprintf(stderr,"too many adjacency entries\n");return 2;}
 uint32_t *deg=calloc(nv,sizeof(*deg)); if(!deg){perror("calloc");return 2;}
 for(uint64_t e=0;e<ne;e++){deg[eu[e]]++;deg[ev[e]]++;}
 uint64_t *off=xmalloc((size_t)(nv+1)*8); off[0]=0;for(uint32_t v=0;v<nv;v++)off[v+1]=off[v]+deg[v];
 uint32_t *adj=xmalloc((size_t)(2*ne)*4);uint64_t *cur=xmalloc((size_t)nv*8);memcpy(cur,off,(size_t)nv*8);
 for(uint64_t e=0;e<ne;e++){adj[cur[eu[e]]++]=ev[e];adj[cur[ev[e]]++]=eu[e];}free(cur);
 uint8_t *col=xmalloc(nv); uint64_t *tabu=xmalloc((size_t)nv*4*sizeof(uint64_t)); memset(tabu,0,(size_t)nv*4*sizeof(uint64_t));
 uint32_t *bad=xmalloc((size_t)nv*4),*pos=xmalloc((size_t)nv*4),*conf=xmalloc((size_t)nv*4);uint32_t nbad=0;
 for(uint32_t v=0;v<nv;v++){col[v]=(uint8_t)pick(4);pos[v]=UINT32_MAX;}
 uint64_t conflicts=0;
 for(uint32_t v=0;v<nv;v++){uint32_t c=0;for(uint64_t j=off[v];j<off[v+1];j++)c+=(col[adj[j]]==col[v]);conf[v]=c;if(c){bad[nbad]=v;pos[v]=nbad++;}}
 for(uint64_t e=0;e<ne;e++)conflicts+=(col[eu[e]]==col[ev[e]]);
 fprintf(stderr,"graph nv=%u ne=%"PRIu64" seed=%"PRIu64" start_conflicts=%"PRIu64"\n",nv,ne,seed,conflicts);
 uint64_t steps=0, since_restart=0, best=conflicts; uint32_t stale=0;
 while(conflicts && steps<maxflips){
  if(nbad==0){fprintf(stderr,"internal bad-list inconsistency\n");return 3;}
  uint32_t v=bad[pick(nbad)], old=col[v], counts[4]={0,0,0,0};
  for(uint64_t j=off[v];j<off[v+1];j++)counts[col[adj[j]]]++;
  uint32_t bestc=UINT32_MAX, choices[4], nc=0;
  for(uint32_t c=0;c<4;c++)if(c!=old && (tabu[4*(size_t)v+c] <= steps || counts[c]==0)){
    if(counts[c]<bestc){bestc=counts[c];nc=0;} if(counts[c]==bestc)choices[nc++]=c;
  }
  if(!nc){for(uint32_t c=0;c<4;c++)if(c!=old){if(counts[c]<bestc){bestc=counts[c];nc=0;}if(counts[c]==bestc)choices[nc++]=c;}}
  uint32_t nw=choices[pick(nc)];
  /* Tabu is stored modulo 256; short tenure keeps the byte representation exact. */
  tabu[4*(size_t)v+old]=steps+5+pick(8);
  col[v]=(uint8_t)nw;
  for(uint64_t j=off[v];j<off[v+1];j++){
   uint32_t u=adj[j], before=conf[u];
   if(col[u]==old) conf[u]--;
   if(col[u]==nw) conf[u]++;
   if(!before&&conf[u]){pos[u]=nbad;bad[nbad++]=u;} else if(before&&!conf[u]){uint32_t p=pos[u],last=bad[--nbad];bad[p]=last;pos[last]=p;pos[u]=UINT32_MAX;}
  }
  uint32_t cv=0;for(uint64_t j=off[v];j<off[v+1];j++)cv+=(col[adj[j]]==nw);conf[v]=cv;
  if(cv==0){uint32_t p=pos[v],last=bad[--nbad];bad[p]=last;pos[last]=p;pos[v]=UINT32_MAX;}
  else { /* v stays conflicted; its list position is unchanged */ }
  /* Update total conflicts from the incident edges whose endpoint is v. */
  conflicts = conflicts - counts[old] + counts[nw];
  steps++; since_restart++;
  if(steps%1000000==0)fprintf(stderr,"step=%"PRIu64" conflicts=%"PRIu64" best=%"PRIu64"\n",steps,conflicts,best);
  if(conflicts<best){best=conflicts;stale=0;}else stale++;
  if(conflicts==0)break;
  if(restartflips&&since_restart>=restartflips){for(uint32_t q=0;q<nv;q++){col[q]=(uint8_t)pick(4);tabu[4*(size_t)q]=tabu[4*(size_t)q+1]=tabu[4*(size_t)q+2]=tabu[4*(size_t)q+3]=0;}
    nbad=0;memset(pos,0xff,(size_t)nv*4);conflicts=0;
    for(uint32_t q=0;q<nv;q++){uint32_t c=0;for(uint64_t j=off[q];j<off[q+1];j++)c+=(col[adj[j]]==col[q]);conf[q]=c;if(c){pos[q]=nbad;bad[nbad++]=q;}conflicts+=c;}
    conflicts/=2;since_restart=0;if(conflicts<best)best=conflicts;fprintf(stderr,"restart step=%"PRIu64" conflicts=%"PRIu64"\n",steps,conflicts);
  }
 }
 if(conflicts==0){for(uint64_t e=0;e<ne;e++)if(col[eu[e]]==col[ev[e]]){fprintf(stderr,"verification failed at edge %"PRIu64"\n",e);return 3;}FILE*g=fopen(out,"wb");if(!g){perror(out);return 2;}if(fwrite(col,1,nv,g)!=nv){fprintf(stderr,"write failed\n");return 2;}fclose(g);printf("SAT seed=%"PRIu64" steps=%"PRIu64" output=%s\n",seed,steps,out);return 0;}
 fprintf(stderr,"UNKNOWN seed=%"PRIu64" steps=%"PRIu64" conflicts=%"PRIu64" best=%"PRIu64"\n",seed,steps,conflicts,best);return 1;
}
