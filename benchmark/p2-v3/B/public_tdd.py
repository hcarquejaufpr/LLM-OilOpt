#!/usr/bin/env python3
import json, subprocess, sys, tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOLVE=ROOT/"solve.py"
PROBLEM=ROOT/"shared"/"problema_p2_public.json"
REQUIRED={"status","epsilon_FWIT","requested","realized","R","FOPT","FWIT","solver_gap"}

def interp(x,spec):
    xs,ys=spec["nodes"],spec["values"]
    assert xs[0]-1e-8 <= x <= xs[-1]+1e-8
    if x <= xs[0]: return ys[0]
    if x >= xs[-1]: return ys[-1]
    for i in range(len(xs)-1):
        if xs[i] <= x <= xs[i+1]:
            t=(x-xs[i])/(xs[i+1]-xs[i])
            return ys[i]+t*(ys[i+1]-ys[i])
    raise AssertionError("PWL interpolation failed")

def main():
    assert PROBLEM.is_file(), "problema público ausente"
    p=json.loads(PROBLEM.read_text())
    assert p["decision_variables"]["fixed"]["INJECT1"] == 79.5
    assert p["decision_variables"]["capacity_requested"] == 636.0
    assert p["normalized_features"]["real_i"] == "(r_i - 80)/20"
    assert len(p["optimization"]["epsilon_candidates"]) == 5
    assert SOLVE.is_file(), "solve.py ausente"

    eps=p["optimization"]["epsilon_candidates"][2]["epsilon_FWIT"]
    with tempfile.TemporaryDirectory() as td:
        out=Path(td)/"solution.json"
        cp=subprocess.run([sys.executable,str(SOLVE),"--epsilon",str(eps),"--output",str(out)],
                          cwd=ROOT,capture_output=True,text=True,timeout=120)
        assert cp.returncode == 0, cp.stderr[-2000:]
        assert out.is_file(), "solve.py não criou o output"
        s=json.loads(out.read_text())
        assert REQUIRED <= set(s), f"campos ausentes: {sorted(REQUIRED-set(s))}"
        assert abs(float(s["epsilon_FWIT"])-eps) <= 1e-6
        assert set(s["requested"]) == set(p["wells"])
        assert set(s["realized"]) == set(p["wells"])
        q={k:float(v) for k,v in s["requested"].items()}
        r={k:float(v) for k,v in s["realized"].items()}

        assert abs(q["INJECT1"]-79.5) <= 1e-6
        assert sum(q.values()) <= 636.0 + 1e-5
        for w,b in p["decision_variables"]["domains"].items():
            assert b[0]-1e-6 <= q[w] <= b[1]+1e-6
        for w in p["wells"]:
            assert r[w] <= q[w] + 1e-6
        R=sum(r.values())
        assert abs(R-float(s["R"])) <= 1e-5
        assert float(s["FWIT"]) <= eps + 1e-3
        # clipping exato: r_i = q_i - clip(deficit_i,0,q_i), deficit_i dos modelos publicos
        W=p["wells"]; Q=(sum(q.values())-600.0)/20.0
        for w in W:
            m=p["layer1_deficit_models"][w]
            if not m: d=0.0
            else:
                own=(q[w]-80.0)/20.0
                d=m.get("intercept",0.0)+m.get("own",0.0)*own+m.get("Q",0.0)*Q
                if w=="INJECT3":
                    d+=m["own2"]*interp(own,p["pwl_functions"]["own2_INJECT3"])+m["Q2"]*interp(Q,p["pwl_functions"]["Q2"])
            assert abs(r[w]-(q[w]-min(max(d,0.0),q[w]))) <= 1e-3, f"clipping inexato em {w}"

        Z=(R-600.0)/20.0
        Z2=interp(Z,p["pwl_functions"]["Z2"])
        fc=p["layer2"]["FOPT"]["coefficients"]
        fcalc=fc["intercept"]
        for i,w in enumerate(p["wells"],1):
            real_i=(r[w]-80.0)/20.0
            fcalc += fc[f"r{i}"]*real_i
        fcalc += fc["Z2"]*Z2
        assert abs(float(s["FOPT"])-fcalc) <= 1.0, (
            f"FOPT incoerente com definição pública v3: {s['FOPT']} vs {fcalc}"
        )
    print("PUBLIC_TDD_V3_PASS")

if __name__=="__main__":
    main()
