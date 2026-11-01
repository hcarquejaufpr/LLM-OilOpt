# Formulação MILP — P2 público v3

Fonte única dos fatos: `shared/problema_p2_public.json`. O programa lê esse arquivo a cada execução. A solução é obtida com `scipy.optimize.milp` (HiGHS).

## Variáveis e camada requested

Para cada poço (i\in\{1,\ldots,8\}), (q_i) usa exatamente o domínio público; (q_1=79.5). Impõe-se

\[
\sum_i q_i\le636,\qquad Q=(\sum_i q_i-600)/20,\qquad own_i=(q_i-80)/20.
\]

O déficit (d_i) é a expressão afim de `layer1_deficit_models`, usando `own2_INJECT3` e `Q2` onde presentes. Modelos vazios dão (d_i=0). Não se impõem limites de `deficit_bounds` sobre (d_i).

## Interpolação pública

Para cada função PWL, com nós (x_j) e valores (y_j), usam-se pesos λ_j e binárias de segmento (s_k):

\[
\sum_j\lambda_j=1,\quad \sum_k s_k=1,\quad
\lambda_j\ge0,\quad
\lambda_j\le s_{j-1}+s_j,\quad
x=\sum_j x_j\lambda_j,\quad y=\sum_j y_j\lambda_j.
\]

Nos extremos, somente o segmento adjacente existe. Os nós e valores são lidos sem alteração do JSON público. A construção vale para `own2_INJECT3`, `Q2` e `Z2`.

## Clipping exato

Para cada poço com modelo de déficit não vazio, três binárias (b_i^-,b_i^0,b_i^+\) somam 1 e ativam, respectivamente:

| Regime | Condições ativas |
|---|---|
| Negativo | (d_i\le0,\ r_i=q_i) |
| Intermediário | (0\le d_i\le q_i,\ r_i=q_i-d_i) |
| Saturado | (d_i\ge q_i,\ r_i=0) |

Cada condição ativa vira uma implicação Big-M. Quando dois regimes se encontram, ambos descrevem o mesmo (r_i). Nos modelos vazios, (d_i=0\) e (r_i=q_i\) diretamente. O limite (0\le r_i\le q_i^{max}\) é uma consequência universal da definição de clipping e serve para limitar as implicações desativadas; `realized_bounds` não são impostos.

Se ([L_i,U_i]) é `deficit_bounds` e ([a_i,b_i]) é o domínio público de (q_i), os M usados são:

| Implicação | M suficiente |
|---|---|
| (d_i\le0\) | \(\max(0,U_i)\) |
| \(r_i-q_i=0\) | \(b_i\) |
| \(d_i\ge0\) | \(\max(0,-L_i)\) |
| \(d_i-q_i\le0\) | \(\max(0,U_i-a_i)\) |
| \(r_i-q_i+d_i=0\) | \(\max(|L_i-b_i|,|b_i-a_i+U_i|)\) |
| \(d_i-q_i\ge0\) | \(\max(0,b_i-L_i)\) |
| \(r_i=0\) | \(b_i\) |

O código soma (10^{-8}\) a cada M para acomodar arredondamentos dos limites auxiliares. Esses limites participam apenas do cálculo de M; não criam restrições de factibilidade para (d_i) ou (r_i).

## Camada realized e objetivo

\[
R=\sum_i r_i,\quad Z=(R-600)/20,\quad real_i=(r_i-80)/20.
\]

`Z2` é a interpolação pública em (Z). FOPT é `intercept + sum(r_i coefficient * real_i) + Z2 coefficient * Z2`, conforme `layer2.FOPT`. FWIT é `intercept + Z coefficient * Z`, conforme `layer2.FWIT`. Maximiza-se FOPT sujeito a `FWIT <= epsilon_FWIT`, com epsilon vindo de `--epsilon`. WBHP não integra o MILP.

## Execução e auditoria

`python3 solve.py --epsilon 1943993.810755 --output resultado.json`

O JSON de saída inclui os campos contratuais, features, déficits, constantes Big-M, gap do solver e erros de recomputação independentes das equações públicas. Status `optimal` só é emitido quando HiGHS retorna ótimo e a auditoria numérica passa.
