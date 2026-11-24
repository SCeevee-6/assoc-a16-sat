import json
import pathlib
import sys

source, target = map(pathlib.Path, sys.argv[1:3])
shard = int(sys.argv[3])
assert shard in range(5)
base = pathlib.Path(__file__).resolve().parent
meta = json.loads((source / "kernel-metadata.json").read_text())
builder = (base.parent / "build_a16.c").read_text()
values = dict(BUILDER=builder, ORACLE=(base/"oracle.cpp").read_text(),
              REFERENCE=(base/"triangle_bits.txt").read_text().strip(), SHARD=shard)
worker = "".join(f"{k} = {v!r}\n" for k,v in values.items()) + (base/"worker.py").read_text()
compile(worker, "run_worker.py", "exec")
(target / "run_worker.py").write_text(worker)
meta.update(code_file="run_worker.py", enable_internet=True)
(target / "kernel-metadata.json").write_text(json.dumps(meta,indent=2))
