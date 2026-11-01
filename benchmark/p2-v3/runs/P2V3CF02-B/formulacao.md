# Formulação MILP — P2 público v3

Fonte única dos dados: `shared/problema_p2_public.json`. O programa `solve.py` lê esse arquivo em cada execução. Requer Python 3, NumPy e SciPy com `scipy.optimize.milp` (HiGHS).

## Variáveis e conjunto factível

Para os oito poços, `q_i` tem os limites de `decision_variables`: `q_1 = 79.5`, e os demais intervalos são lidos literalmente do JSON. Impõe-se `Σ_i q_i ≤ 636`. As variáveis `r_i` têm apenas os limites lógicos `0 ≤ r_i ≤ max(q_i)`; `r_i ≤ q_i` também é consequência do clipping. O déficit `d_i` é uma variável livre vinculada às equações públicas, sem limite próprio. Nenhum limite de `deficit_bounds` ou `realized_bounds` é imposto a `d_i` ou `r_i`.

As features são equações lineares exatas:

```
20 Q       = Σ_i q_i − 600
20 own_i   = q_i − 80
R          = Σ_i r_i
20 Z       = R − 600
20 real_i  = r_i − 80
```

O déficit de cada poço é `intercept + βown·own_i + βQ·Q + βown2·own2_INJECT3 + βQ2·Q2`, usando somente os termos presentes em `layer1_deficit_models`. Para um modelo vazio, `d_i = 0` e `r_i = q_i`.

## Interpolação por segmentos adjacentes

Para cada uma das três funções `own2_INJECT3`, `Q2` e `Z2`, sejam `(x_j,v_j)` os nós e valores **literais** de `pwl_functions`. Para cada segmento `j`, usam-se `y_j ∈ {0,1}` e `0 ≤ t_j ≤ y_j`, com:

```
Σ_j y_j = 1
x = Σ_j [x_j y_j + (x_{j+1} − x_j)t_j]
f = Σ_j [v_j y_j + (v_{j+1} − v_j)t_j]
```

Assim exatamente um segmento adjacente é escolhido e a interpolação linear pública é reproduzida. Os próprios nós definem o domínio da função interpolada.

## Clipping exato em três regimes

Para cada poço com modelo de déficit, `z0+z1+z2=1`, `z_k∈{0,1}`. Os regimes são:

| Regime | Condição | Vazão realizada |
| --- | --- | --- |
| `z0=1` | `d≤0` | `r=q` |
| `z1=1` | `0≤d≤q` | `r=q−d` |
| `z2=1` | `d≥q` | `r=0` |

Sejam `[L,U]` os `deficit_bounds` auxiliares e `[l,u]` o domínio de `q`. Todos os `M` abaixo são coeficientes de relaxação das implicações; **não** aparecem como limites de factibilidade sobre `d` ou `r`. Define-se `A=max(0,U)`, `B=max(0,−L)`, `C=max(0,U−l)` e `D=max(0,u−L)`. Com `0≤r≤q`, as desigualdades implementadas são:

```
d + A z0                  ≤ A
−r + q + u z0             ≤ u
−d + B z1                 ≤ B
d − q + C z1              ≤ C
r − q + d + A z1          ≤ A
−r + q − d + D z1         ≤ D
−d + q + D z2             ≤ D
r + u z2                  ≤ u
```

As constantes cobrem todas as implicações desligadas: `d≤U`, `d≥L`, `d−q≤U−l`, `d−q≥L−u`, `r−q∈[−u,0]` e `r−q+d∈[L−u,U]`. Os bounds auxiliares são usados somente nessas deduções. Nos pontos de empate (`d=0` ou `d=q`), mais de um regime pode ser válido, mas todos geram o mesmo `r`.

## Camada 2 e objetivo

`FOPT = intercept + Σ_i β_i real_i + β_Z2 Z2`, com os coeficientes de `layer2.FOPT`; `FWIT = intercept + β_Z Z`, com os coeficientes de `layer2.FWIT`. Maximiza-se `FOPT` sujeito a `FWIT≤epsilon_FWIT`, onde `epsilon_FWIT` vem de `--epsilon`. WBHP é validação física externa e não integra o MILP. Nenhuma solução de referência ou resultado OPM é usado.

## Execução

```
python3 solve.py --epsilon 1943993.810755 --output resultado.json
python3 B/public_tdd.py
```

O resultado inclui vazões solicitadas e realizadas, déficits, features, métricas, constantes Big-M e estado/gap do HiGHS para auditoria.
