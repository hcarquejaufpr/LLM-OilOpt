# Condição B — SDD do benchmark P2 v3

Implemente o mesmo problema descrito em `shared/problema_p2_public.json`, seguindo os requisitos abaixo.

## Requisitos

- **R01** — Representar `q_i` para INJECT1–INJECT8, com INJECT1 fixo e demais domínios exatamente como especificados.
- **R02** — Impor `sum(q_i) <= 636` sobre vazões solicitadas.
- **R03** — Calcular exatamente as features públicas: `Q=(sum(q)-600)/20`, `own_i=(q_i-80)/20`, `real_i=(r_i-80)/20` e `Z=(sum(r)-600)/20`.
- **R04** — Construir os modelos de déficit da camada requested -> realized usando somente os coeficientes públicos.
- **R05** — Representar `r_i = q_i - clip(deficit_i,0,q_i)` exatamente por três regimes. Se usar Big-M, derivar constantes suficientes para todas as implicações a partir de `deficit_bounds` e dos domínios de `q`.
- **R06** — Implementar `own2_INJECT3`, `Q2` e `Z2` por SOS2 ou MILP equivalente com adjacência, preservando os nós públicos.
- **R07** — Definir `R=sum(r_i)` e `Z=(R-600)/20`.
- **R08** — Calcular FOPT usando os coeficientes públicos sobre `real_i=(r_i-80)/20` e o termo público `Z2`; calcular FWIT exatamente segundo a camada 2 pública.
- **R09** — Maximizar FOPT com `FWIT <= epsilon_FWIT`, sendo epsilon parâmetro de execução.
- **R15** — `deficit_bounds` are auxiliary bounds supplied solely for deriving and validating valid Big-M constants for the exact clipping formulation. They MUST NOT be imposed as additional feasibility constraints on `deficit_i`. `realized_bounds` are auxiliary bounds on `r_i`, supplied solely for deriving and validating Big-M constants. They MUST NOT be imposed as additional feasibility constraints on `r_i`. O conjunto factível é exatamente o descrito em `feasible_set`.
- **R10** — Não incluir WBHP como restrição do MILP.
- **R11** — Entregar `formulacao.md`, `solve.py`, `model.json`, `evidence.json`.
- **R12** — `solve.py` deve aceitar `--epsilon <valor> --output <arquivo.json>` e produzir resultado auditável.
- **R13** — `evidence.json` deve mapear R01–R12 para implementação e testes/evidências.
- **R14** — Não usar solução de referência ou resultados OPM.

## Critério de conclusão

A entrega está concluída quando R01–R15 estiverem implementados e os testes públicos aplicáveis estiverem verdes. Não altere os fatos do problema para fazer os testes passarem.
