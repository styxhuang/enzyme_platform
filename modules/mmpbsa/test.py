import os
import sys
import json

def _import_func():
    try:
        from app import _parse_final_decomp_dat
        return _parse_final_decomp_dat
    except Exception:
        from modules.mmpbsa.app import _parse_final_decomp_dat
        return _parse_final_decomp_dat

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error":"usage: python test.py <FINAL_DECOMP_MMPBSA.dat>"}))
        sys.exit(1)
    inp = sys.argv[1]
    if not os.path.exists(inp):
        print(json.dumps({"error":f"not found: {inp}"}))
        sys.exit(1)
    outdir = os.path.dirname(os.path.abspath(inp))
    func = _import_func()
    res = {}
    try:
        res = func(inp, outdir)
    except Exception as e:
        print(json.dumps({"error":str(e)}))
        sys.exit(1)
    files = []
    try:
        for n in os.listdir(outdir):
            p = os.path.join(outdir, n)
            if os.path.isfile(p):
                files.append(n)
    except Exception:
        pass
    print(json.dumps({"outputs":res, "files":files}, ensure_ascii=False))

if __name__ == "__main__":
    main()

