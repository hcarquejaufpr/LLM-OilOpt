# Formulação MILP — P2 público v3

Fonte única dos parâmetros: `shared/problema_p2_public.json`. O programa lê o arquivo em cada execução. Não usa resultados de referência nem OPM.

## Variáveis e equações

Para cada poço `i`, `q_i` tem o domínio público e `q_1=79.5`. Impõe-se `sum(q_i)<=636`. As features são igualdades:

```text
20 Q = sum(q_i)-600
20 own_i = q_i-80
20 real_i = r_i-80
R = sum(r_i)
20 Z = R-600
```

O déficit `d_i` é a combinação linear dos coeficientes públicos de `layer1_deficit_models`. Coeficientes ausentes são zero. Para INJECT3, os termos `own2_INJECT3` e `Q2` são as funções lineares por partes públicas. Déficits vazios de INJECT2, INJECT4 e INJECT5 equivalem a zero.

Cada função PWL (`own2_INJECT3`, `Q2`, `Z2`) usa todos os nós e valores públicos. Para os nós `x_j,y_j`, existem `lambda_j>=0`, `sum(lambda)=1`, e binárias de segmento `s_k` com `sum(s)=1`. As equações são `x=sum(x_j lambda_j)`, `y=sum(y_j lambda_j)`. Em cada nó, `lambda_j <= s_{j-1}+s_j`, omitindo segmentos inexistentes nas pontas. Assim somente dois nós adjacentes podem ter peso positivo. Os intervalos das funções PWL coincidem com os limites implícitos do conjunto factível público.

Para cada poço, três binárias `a,b,c` somam 1. Elas selecionam exatamente um regime:

| Regime | Condições ativas |
| --- | --- |
| `a` | `d<=0`, `r=q` |
| `b` | `0<=d<=q`, `r=q-d` |
| `c` | `d>=q`, `r=0` |

Esses regimes dão `r=q-clip(d,0,q)`, inclusive nas fronteiras. `d` e `r` são variáveis sem limites auxiliares de domínio: `deficit_bounds` e `realized_bounds` não são restrições de factibilidade.

## Constantes Big-M

Sejam `[L_d,U_d]` os `deficit_bounds` auxiliares públicos e `[L_q,U_q]` o domínio público de `q`. As constantes abaixo relaxam somente as implicações dos regimes não selecionados:

```text
a => d<=0:       M=max(0,U_d)
b => d>=0:       M=max(0,-L_d)
b => d-q<=0:     M=max(0,U_d-L_q)
c => d-q>=0:     M=max(0,U_q-L_d)
a => r-q=0:      M=U_q em ambas as direções
b => r-q+d=0:    M_positivo=max(0,U_d-L_q), M_negativo=max(0,-L_d)
c => r=0:        M=U_q em ambas as direções
```

Por exemplo, `d<=M(1-a)`, `d>=-M(1-b)`, `r-q+d<=M_positivo(1-b)` e `r-q+d>=-M_negativo(1-b)`. Os limites para as igualdades seguem da união dos três regimes: `r-q` está em `[-U_q,0]`, `r-q+d` em `[L_d,max(0,U_d-L_q)]`, e `r` em `[0,U_q]`. Esses limites são usados para provar validade das constantes; não são adicionados como limites de variável. `realized_bounds` também são apenas auxiliares e não são impostos.

## Objetivo

```text
FOPT = intercept_FOPT + sum_i beta_i real_i + beta_Z2 Z2
FWIT = intercept_FWIT + beta_Z Z
maximize FOPT
subject to FWIT <= epsilon_FWIT
```

`epsilon_FWIT` vem de `--epsilon`. WBHP pertence somente à validação física externa e não entra no MILP.

O programa usa `scipy.optimize.milp` (HiGHS). O resultado inclui vazões solicitadas e realizadas, déficits recalculados, regimes, FOPT, FWIT, gap do solver e resíduos de auditoria. A auditoria recalcula clipping, FOPT e FWIT diretamente do JSON público e falha se houver inconsistência.
