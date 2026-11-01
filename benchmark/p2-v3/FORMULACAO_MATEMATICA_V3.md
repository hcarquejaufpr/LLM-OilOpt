# Formulação matemática v3 (especificação de auditoria)

Variáveis: q_1=79,5; q_i∈dom_i (i=2..8); Σq ≤ 636.
Features: Q=(Σq−600)/20; own_i=(q_i−80)/20.
Déficit: d_i = a_i + b_i·own_i + c_i·Q (+ para INJECT3: e·PWL_own2(own_3) + f·PWL_Q2(Q)); INJECT2/4/5: d=0.
Clipping: r_i = q_i − clip(d_i,0,q_i), exato em 3 regimes (d≤0: r=q; 0<d<q: r=q−d; d≥q: r=0), via binárias e Big-M derivado de `deficit_bounds` e dos domínios de q.
R=Σr_i; Z=(R−600)/20; real_i=(r_i−80)/20.
FOPT = α + Σ β_i·real_i + β_Z2·PWL_Z2(Z). FWIT = γ + δ·Z.
Problema: max FOPT s.a. FWIT ≤ ε (ε ∈ {P1..P5}), Σq ≤ 636, domínios.
PWL por SOS2 ou adjacência equivalente. WBHP externo. `deficit_bounds` e `realized_bounds`: NÃO são restrições.
Fonte numérica: `shared/problema_p2_public.json` (gerado por `GERAR_PROBLEMA_PUBLICO_V3.py` a partir do v2 somente leitura).
