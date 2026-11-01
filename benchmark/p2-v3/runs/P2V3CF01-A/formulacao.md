# Formulação MILP/PWL — P2 público v3

## Dados, decisão e objetivo

A única fonte de coeficientes, domínios e nós é
`shared/problema_p2_public.json` (`p2-public-problem-v3-semantic-closed`).
Para cada poço (i=1,\ldots,8), (q_i) é a vazão solicitada, (d_i) o
déficit calculado, (t_i=\operatorname{clip}(d_i,0,q_i)) e (r_i=q_i-t_i)
a vazão realizada. `R` na saída é (\sum_i r_i). As unidades e os valores
numéricos são os do arquivo público.

Maximiza-se `FOPT` com (q_1=79{,}5), os domínios públicos de
(q_2,\ldots,q_8), (\sum_i q_i\le636) e

\[
\mathrm{FWIT}\le\epsilon_{\mathrm{FWIT}}.
\]

O parâmetro epsilon é lido na execução. Nenhum valor de solução é codificado.

## Features e déficit

As igualdades são

\[
Q=(\textstyle\sum_i q_i-600)/20,\quad
\operatorname{own}_i=(q_i-80)/20,\quad
Z=(\textstyle\sum_i r_i-600)/20,\quad
\operatorname{real}_i=(r_i-80)/20.
\]

Para cada poço, o déficit é a soma do intercepto e dos termos presentes em
`layer1_deficit_models`:

\[
d_i=a_i+b_i^{\rm own}\operatorname{own}_i+b_i^Q Q
   +b_i^{\rm own2}P_{\rm own2}(\operatorname{own}_3)
   +b_i^{Q2}P_{Q2}(Q).
\]

Coeficientes ausentes valem zero. Assim (d_2=d_4=d_5=0). Cada (P)
interpola linearmente os pares `nodes`/`values` públicos, sem extrapolação.

## PWL com adjacência exata

Para cada função (P) com nós ((x_j,y_j)), são usadas variáveis
\(\lambda_j\ge0\) e binárias \(s_k\) para os segmentos entre nós
consecutivos:

\[
\sum_j\lambda_j=1,\quad \sum_k s_k=1,\quad
x=\sum_j x_j\lambda_j,\quad P(x)=\sum_j y_j\lambda_j,
\quad \lambda_j\le s_{j-1}+s_j.
\]

Índices de segmento inexistentes têm valor zero. Logo, apenas os dois
vértices de um segmento escolhido podem ter peso positivo. A construção é
equivalente à adjacência SOS2 e é usada para `own2_INJECT3`, `Q2` e `Z2`.
Os domínios dos nós cobrem os argumentos gerados pelo conjunto factível
público. Os limites de `realized_bounds` servem apenas para conferir essa
cobertura; não aparecem nas variáveis ou restrições do MILP.

## Clipping em três regimes

Para cada poço, (z_i^-,z_i^0,z_i^+\in\{0,1\}) e
(z_i^-+z_i^0+z_i^+=1). Os regimes selecionados impõem exatamente:

| Regime | Condição | Déficit cortado |
| --- | --- | --- |
| \(z_i^-=1\) | \(d_i\le0\) | \(t_i=0\) |
| \(z_i^0=1\) | \(0\le d_i\le q_i\) | \(t_i=d_i\) |
| \(z_i^+=1\) | \(d_i\ge q_i\) | \(t_i=q_i\) |

Depois, (r_i=q_i-t_i). Nos pontos de encontro, qualquer regime vizinho
produz o mesmo (r_i). O limite interno (0\le t_i\le\bar q_i) decorre
da definição do clipping. `realized_bounds` não é usado como limite de
(r_i), que é uma variável contínua sem limites próprios.

Os Big-M são derivados por poço de
(L_i=\underline d_i-10^{-7}), (U_i=\overline d_i+10^{-7}) e
(\underline q_i,\bar q_i). O pequeno alargamento cobre arredondamento
decimal. As implicações implementadas são:

\[
\begin{aligned}
d_i &\le U_i^+(1-z_i^-), &t_i&\le\bar q_i(1-z_i^-),\\
d_i &\ge L_i^-(1-z_i^0),
&d_i-q_i&\le (U_i-\underline q_i)^+(1-z_i^0),\\
t_i-d_i&\le(\bar q_i-L_i)(1-z_i^0),
&d_i-t_i&\le U_i^+(1-z_i^0),\\
d_i-q_i&\ge-(\bar q_i-L_i)^+(1-z_i^+),
&t_i-q_i&\le(\bar q_i-\underline q_i)(1-z_i^+),\\
q_i-t_i&\le\bar q_i(1-z_i^+).
\end{aligned}
\]

Aqui (a^+=\max(a,0)) e (L_i^-=\min(L_i,0)). Essas constantes
cobrem todos os valores inativos possíveis, dados os limites auxiliares
válidos de déficit e os domínios de (q_i). **`deficit_bounds` só entra no
cálculo dos Big-M; nenhuma restrição (L_i\le d_i\le U_i) é adicionada.**
Os valores calculados por poço constam de `model.json`.

## Camada final

Com (P_{Z2}(Z)) construído dos nós públicos,

\[
\begin{aligned}
\mathrm{FOPT}&=501555.5938340376
 +\sum_{i=1}^{8}\beta_i\operatorname{real}_i
 -167.73716538798362P_{Z2}(Z),\\
\mathrm{FWIT}&=2159936.0100392425+71980.73309472566Z.
\end{aligned}
\]

Os oito \(\beta_i\) são lidos de `layer2.FOPT.coefficients`; não são
substituídos por coeficientes sobre (r_i). A pressão WBHP e seu limite
de 420 bar pertencem exclusivamente à validação física externa e não
participam do MILP.

## Execução e auditoria

`python3 solve.py --epsilon <valor> --output <arquivo.json>` usa
`scipy.optimize.milp`/HiGHS. A saída registra o estado, as vazões, `R`,
FOPT, FWIT, gap, déficits e resíduos de uma reavaliação independente das
fórmulas públicas. `evidence.json` reúne os cinco epsilons públicos.
