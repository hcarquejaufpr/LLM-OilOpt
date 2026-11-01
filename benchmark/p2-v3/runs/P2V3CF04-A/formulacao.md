# Formulação MILP/PWL — P2 público v3

## Dados e objetivo

Os únicos dados numéricos do modelo vêm de `shared/problema_p2_public.json`.
Para cada poço (i=1,\ldots,8), (q_i) é a vazão solicitada, (d_i) é o
déficit calculado pela primeira camada e (r_i) é a vazão realizada. Maximiza-se
`FOPT` sujeito a `FWIT <= epsilon_FWIT`. O horizonte de 3600 dias já está
incorporado nos coeficientes públicos; não se multiplica a função objetivo
ou a restrição por esse número.

O conjunto de decisão para as vazões solicitadas é (q_1=79.5), os sete
intervalos públicos de `decision_variables.domains` e
\(\sum_i q_i\le636\). Não se impõe WBHP no MILP.

## Features e déficits

\[
Q=(\sum_i q_i-600)/20,\quad own_i=(q_i-80)/20,\\
Z=(\sum_i r_i-600)/20,\quad real_i=(r_i-80)/20.
\]

Para cada poço, usando coeficientes ausentes como zero:

\[
d_i=a_i+b_i^Q Q+b_i^{own}own_i
  +b_i^{own2}P_{own2}(own_3)+b_i^{Q2}P_{Q2}(Q).
\]

Os dois últimos termos aparecem somente em `INJECT3`. Cada (P) interpola
linearmente os nós e valores da função pública; em particular, não se avalia
um quadrado contínuo entre nós.

## PWL com adjacência exata

Para cada função PWL com nós ordenados \((x_j,y_j)\), introduzem-se
\(\lambda_j\ge0\) e variáveis binárias \(s_j\) para os segmentos:

\[
\sum_j\lambda_j=1,\quad\sum_j s_j=1,\quad
x=\sum_j x_j\lambda_j,\quad y=\sum_j y_j\lambda_j,\quad
\lambda_j\le s_{j-1}+s_j,
\]

omitindo os segmentos inexistentes nas extremidades. Portanto somente
os dois nós de um mesmo segmento podem receber peso. Aplica-se a
`own2_INJECT3`, `Q2` e `Z2`.

## Clipping em três regimes

Defina \(h_i=\max(d_i,0)\), \(g_i=\max(d_i-q_i,0)\) e
\(r_i=q_i-h_i+g_i\). Como (q_i>0), isto equivale exatamente a:

| Regime | Condição | Vazão realizada |
|---|---|---|
| 1 | \(d_i\le0\) | \(r_i=q_i\) |
| 2 | \(0\le d_i\le q_i\) | \(r_i=q_i-d_i\) |
| 3 | \(d_i\ge q_i\) | \(r_i=0\) |

Logo \(r_i=q_i-\operatorname{clip}(d_i,0,q_i)\), incluindo as fronteiras.
Os dois operadores `max` são linearizados exatamente. Para um argumento
\(t\in[L,U]\), (v=\max(t,0)\), binária (b), usam-se:

\[
v\ge t,\quad v\ge0,\quad
v\le t-L(1-b),\quad v\le Ub.
\]

No primeiro operador, \(t=d_i\), \(L=D_i^-\), \(U=D_i^+\) vêm de
`deficit_bounds`. No segundo, \(t=d_i-q_i\), usam-se os limites seguros
\(L=D_i^- - q_i^{max}\) e \(U=D_i^+ - q_i^{min}\). Estes Big-M são válidos
para todas as implicações, inclusive quando (U<0), situação em que a
binária correspondente fica forçada ao regime inativo. Os limites de
`deficit_bounds` são usados **somente como constantes** dessas desigualdades:
o modelo não contém limites de variável nem restrições diretas
\(D_i^-\le d_i\le D_i^+\). Similarmente, não se impõe `realized_bounds`
a (r_i); estes limites não são necessários para esta linearização.
Nenhuma restrição de WBHP é criada.

## Segunda camada

Com os coeficientes públicos `layer2`:

\[
FOPT=a_F+\sum_{i=1}^{8}\beta_i^{F}real_i
       +\beta_{Z2}^{F}P_{Z2}(Z),\qquad
FWIT=a_W+\beta_Z^W Z.
\]

O objetivo é maximizar (FOPT\), sujeito a (FWIT\le\epsilon_{FWIT}\).
`R` na saída é \(\sum_i r_i\). O mesmo modelo é reconstruído a partir do
arquivo público para qualquer epsilon finito informado na linha de comando.

## Implementação e tolerâncias

`solve.py` monta o MILP com `scipy.optimize.milp`/HiGHS. O gap relativo
solicitado é \(10^{-9}\); `solver_gap` é o gap reportado pelo solver. A
interpolação PWL e a integralidade são respeitadas dentro das tolerâncias
numéricas do HiGHS. `evidence.json` registra resultados dos cinco epsilon
públicos e resíduos recalculados independentemente.
