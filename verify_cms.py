import json
import struct
import numpy as np

with open("A16_graph.bin", "rb") as f:
    nv = struct.unpack("<i", f.read(4))[0]
    ne = struct.unpack("<q", f.read(8))[0]
    codes = np.frombuffer(f.read(nv * 4), dtype=np.uint32)
    u = np.frombuffer(f.read(ne * 4), dtype=np.int32)
    v = np.frombuffer(f.read(ne * 4), dtype=np.int32)

status = "RUNNING_OR_UNKNOWN"
pos_lits = []
with open("kissat_a16.log", "r") as f:
    for line in f:
        if line.startswith("s "):
            status = line.strip()
        elif line.startswith("v "):
            for tok in line.split()[1:]:
                val = int(tok)
                if val > 0:
                    pos_lits.append(val - 1)

print("Current status:", status)

verified = False
if "SATISFIABLE" in status and "UNSATISFIABLE" not in status:
    pos = np.array(pos_lits, dtype=np.int32)
    colors = np.full(nv, -1, dtype=np.int8)
    colors[pos // 4] = (pos % 4).astype(np.int8)
    verified = bool(np.all((colors >= 0) & (colors <= 3)) and np.all(colors[u] != colors[v]))
    if verified:
        np.savez_compressed("A16_4coloring_verified.npz", codes=codes, colors=colors)
    if verified:
        np.savez_compressed("A16_4coloring_verified.npz", codes=codes, colors=colors)
        print("SUCCESS: CHI(A_16) = 4 VERIFIED!")
    else:
        print("ERROR: Failed edge check!")
elif "UNSATISFIABLE" in status:
    print("SUCCESS: CHI(A_16) >= 5 PROVEN!")

with open("A16_result.json", "w") as f:
    json.dump({"n": 16, "V": int(nv), "E": int(ne), "status": status, "verified": verified}, f, indent=2)
