#!/usr/bin/env python3
"""Compile and launch reproducible direct-SLS workers for A16_graph.bin.

Works in Kaggle notebooks and Codespaces. It does not call the Kaggle API;
start one copy in each prepared Kaggle kernel after making these files and the
graph available in that kernel. A found coloring stops sibling workers.
"""
import argparse, json, os, pathlib, signal, subprocess, sys, time

HERE = pathlib.Path(__file__).resolve().parent

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--graph', required=True, help='path to A16_graph.bin')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--kernel-id', type=int, default=0, help='0..4; reserve a disjoint seed block per Kaggle kernel')
    ap.add_argument('--start-seed', type=int, help='override kernel-based seed block')
    ap.add_argument('--max-flips', type=int, default=100_000_000)
    ap.add_argument('--restart-flips', type=int, default=10_000_000)
    ap.add_argument('--out-dir', default='batch2_out')
    ap.add_argument('--compiler', default=os.environ.get('CC', 'cc'))
    a = ap.parse_args()
    if not (1 <= a.workers <= 20): ap.error('--workers must be 1..20')
    if not (0 <= a.kernel_id <= 99): ap.error('--kernel-id must be 0..99')
    graph=pathlib.Path(a.graph).resolve()
    if not graph.is_file(): ap.error(f'graph not found: {graph}')
    out=pathlib.Path(a.out_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    exe=out/'a16_sls'
    subprocess.run([a.compiler,'-O3','-DNDEBUG','-std=c11','-Wall','-Wextra',str(HERE/'a16_sls.c'),'-o',str(exe)],check=True)
    base=a.start_seed if a.start_seed is not None else 202611230000 + a.kernel_id*1000
    manifest={'graph':str(graph),'workers':a.workers,'kernel_id':a.kernel_id,'seed_base':base,
              'max_flips':a.max_flips,'restart_flips':a.restart_flips,'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'runs':[]}
    procs=[]
    for w in range(a.workers):
        seed=base+w
        dest=out/f'coloring_seed_{seed}.colors'
        log=out/f'worker_{w}.log'
        lf=open(log,'wb')
        cmd=[str(exe),str(graph),'--seed',str(seed),'--max-flips',str(a.max_flips),'--restart-flips',str(a.restart_flips),'--out',str(dest)]
        proc=subprocess.Popen(cmd,stdout=lf,stderr=subprocess.STDOUT)
        procs.append((proc,lf,w,seed,dest,log))
        manifest['runs'].append({'worker':w,'seed':seed,'pid':proc.pid,'output':str(dest),'log':str(log),'status':'RUNNING'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    pending={x[0].pid:x for x in procs}; found=False
    try:
        while pending:
            for pid,item in list(pending.items()):
                proc,lf,w,seed,dest,log=item; code=proc.poll()
                if code is None: continue
                lf.close(); pending.pop(pid)
                manifest['runs'][w]['status']='SAT' if code==0 else ('UNKNOWN' if code==1 else f'ERROR_{code}')
                manifest['runs'][w]['exit_code']=code
                if code==0:
                    found=True; manifest['solution']=str(dest)
                    for other,_,_,_,_,_ in pending.values(): other.terminate()
            (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            if pending: time.sleep(2)
    except KeyboardInterrupt:
        for proc,_,_,_,_,_ in pending.values(): proc.send_signal(signal.SIGINT)
        for proc,lf,_,_,_,_ in pending.values(): proc.wait(); lf.close()
        raise
    manifest['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
    return 0 if found else 1

if __name__=='__main__': sys.exit(main())
