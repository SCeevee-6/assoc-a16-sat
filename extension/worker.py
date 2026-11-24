import itertools
import json
import os
SHARD = int(os.environ.get("A16_SHARD", SHARD))
import pathlib
import subprocess
import sys
import time

# The launcher embeds BUILDER, ORACLE, REFERENCE and SHARD above this file.
out = pathlib.Path(os.environ.get("A16_OUT", "/kaggle/working"))
out.mkdir(parents=True, exist_ok=True)
os.chdir(out)
try:
    import z3
except ImportError:
    subprocess.run([sys.executable, "-m", "pip", "install", "z3-solver==4.13.3.0"], check=True)
    import z3

def report(status, **fields):
    path = pathlib.Path("extension_status.tmp")
    path.write_text(json.dumps(dict(status=status, shard=SHARD, **fields)) + "\n")
    path.replace("extension_status.json")
    print(status, fields, flush=True)

def start(n):
    subprocess.run(["/tmp/ext-builder", str(n), "--graph-only"], check=True)
    return subprocess.Popen(["/tmp/ext-oracle", "A16_graph.bin"], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, text=True)

def ask(proc, bits):
    proc.stdin.write(bits + "\n")
    proc.stdin.flush()
    line = proc.stdout.readline().strip()
    if not line:
        raise RuntimeError("oracle stopped")
    return line

def close(proc):
    proc.stdin.close()
    if proc.wait() != 0:
        raise RuntimeError("oracle failed")

def run():
    pathlib.Path("/tmp/ext-builder.c").write_text(BUILDER)
    pathlib.Path("/tmp/ext-oracle.cpp").write_text(ORACLE)
    subprocess.run(["gcc","-O3","-Werror=int-conversion","/tmp/ext-builder.c","-o","/tmp/ext-builder"],check=True)
    subprocess.run(["g++","-O3","-std=c++17","/tmp/ext-oracle.cpp","-o","/tmp/ext-oracle"],check=True)
    old = dict(zip(itertools.combinations(range(1,15),3),map(int,REFERENCE)))
    assert len(old)==364 and len(REFERENCE)==364
    proc=start(15)
    seed="".join(str(old.get(t,0)) for t in itertools.combinations(range(15),3))
    if ask(proc,seed)!="FOUND": raise RuntimeError("A15 published control failed")
    close(proc)
    pathlib.Path("candidate.colors").unlink()
    proc=start(16)
    # The open oracle has loaded the graph. The output archive needs no graph copy.
    triangles=list(itertools.combinations(range(16),3))
    weights=[z3.BoolVal(False) if t[0]==0 else z3.BoolVal(bool(old[t])) if t in old
             else z3.Bool("t_%d_%d_%d"%t) for t in triangles]
    free=[w for w in weights if not (z3.is_true(w) or z3.is_false(w))]
    assert len(free)==91
    solver=z3.Solver()
    for i in range(5): solver.add(free[i]==bool(SHARD & (1<<i)))
    start_time=time.monotonic()
    rounds=0
    try:
        while time.monotonic()-start_time<36000:
            solver.set(timeout=120000)
            result=solver.check()
            if result!=z3.sat:
                report("SHARD_UNSAT_REQUIRES_REPLAY" if result==z3.unsat else "UNKNOWN", rounds=rounds)
                return
            m=solver.model()
            seed="".join("1" if z3.is_true(m.eval(w,model_completion=True)) else "0" for w in weights)
            line=ask(proc,seed)
            if line=="FOUND":
                pathlib.Path("triangle_seed_16.txt").write_text(seed+"\n")
                report("CANDIDATE_REQUIRES_INDEPENDENT_VERIFICATION",rounds=rounds)
                return
            tag,*raw=line.split()
            if tag!="CYCLE" or len(raw)%4: raise RuntimeError("invalid oracle output")
            indices=list(map(int,raw))
            terms=[]
            for j in range(0,len(indices),4):
                a,b,c,d=(weights[k] for k in indices[j:j+4])
                terms.append(z3.Xor(z3.Xor(a,b),z3.Xor(c,d)))
            solver.add(z3.Or(*terms))
            with open("learned_cycles.jsonl","a") as f: f.write(json.dumps(indices)+"\n")
            rounds+=1
            if rounds==1 or rounds%25==0: report("SEARCHING",rounds=rounds)
        report("TIME_LIMIT",rounds=rounds)
    finally:
        close(proc)
        pathlib.Path("A16_graph.bin").unlink(missing_ok=True)
def child(shard, outdir, logpath):
    env = os.environ.copy()
    env["A16_SHARD"] = str(shard)
    env["A16_OUT"] = str(outdir)
    with open(logpath, "w") as log:
        return subprocess.Popen([sys.executable, __file__], env=env, stdout=log, stderr=subprocess.STDOUT)

if __name__ == "__main__":
    if os.environ.get("A16_CHILD") == "1":
        try:
            run()
        except Exception as exc:
            report("ERROR", error=str(exc))
            raise
    else:
        base = SHARD * 4
        children = []
        for offset in range(4):
            shard = base + offset
            child_out = out / ("shard_%02d" % shard)
            child_out.mkdir(parents=True, exist_ok=True)
            env = os.environ.copy()
            env.update(A16_CHILD="1", A16_SHARD=str(shard), A16_OUT=str(child_out))
            log = open(out / ("shard_%02d.log" % shard), "w")
            proc = subprocess.Popen([sys.executable, __file__], env=env, stdout=log, stderr=subprocess.STDOUT)
            children.append((shard, proc, log))
        results = []
        for shard, proc, log in children:
            code = proc.wait()
            log.close()
            results.append({"shard": shard, "exit_code": code})
        (out / "parallel_status.json").write_text(json.dumps(results, indent=2) + "\\n")
