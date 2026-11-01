# Formulação MILP — P2 público v3

Fonte única dos dados: `shared/problema_p2_public.json`. As variáveis são contínuas, salvo as binárias dos regimes de clipping e dos segmentos PWL. `solve.py` usa `scipy.optimize.milp` (HiGHS). `model.json` contém a matriz linear gerada, com o limite de FWIT marcado como parâmetro.

## Vazões e features

Para cada poço (i), (q_i) tem o domínio público; (q_1=79.5). Impõe-se (sum_iq_i\le636). As igualdades são

\[
Q=(\sum_iq_i-600)/20,\quad own_i=(q_i-80)/20,\quad
R=\sum_i r_i,\quad Z=(R-600)/20,\quad real_i=(r_i-80)/20.
\]

Cada déficit (d_i) é a soma de intercepto e termos presentes em `layer1_deficit_models[i]`, usando `own_i`, `Q`, `own2_INJECT3` e `Q2` conforme os coeficientes públicos. Modelos vazios dão (d_i=0).

## Interpolação com adjacência

Para cada função PWL (f) (`own2_INJECT3`, `Q2`, `Z2`), todos os nós (x_j) e valores (y_j) são lidos sem alteração do JSON público. Com pesos (\lambda_j\ge0) e binárias de segmento (s_k):

\[
\sum_j\lambda_j=1,\quad \sum_ks_k=1,\quad
x=\sum_jx_j\lambda_j,\quad f(x)=\sum_jy_j\lambda_j,\quad
\lambda_j\le s_{j-1}+s_j,
\]

com índices de segmento inexistentes omitidos. Assim, só os dois pesos do segmento escolhido podem ser positivos. A cobertura dos nós de `Q2` e `Z2` decorre dos domínios públicos e dos limites auxiliares válidos, sem introduzir estes últimos como restrições.

## Clipping exato

Para cada poço, binárias (a_i,b_i,c_i) somam 1 e indicam, respectivamente:

| Binária ativa | Condições | Vazão realizada |
|---|---|---|
| (a_i) | (d_i\le0) | (r_i=q_i) |
| (b_i) | (0\le d_i\le q_i) | (r_i=q_i-d_i) |
| (c_i) | (d_i\ge q_i) | (r_i=0) |

As implicações usam (L_i,U_i) de `deficit_bounds` e (l_i,u_i) do domínio de (q_i). Os Big-M suficientes, calculados em `clipping()`, são:

| Implicação | Big-M |
|---|---|
| (a_i\Rightarrow d_i\le0) | \(\max(0,U_i)\) |
| (a_i\Rightarrow r_i=q_i) | \(u_i\) em cada sentido |
| (b_i\Rightarrow d_i\ge0) | \(\max(0,-L_i)\) |
| (b_i\Rightarrow d_i\le q_i) | \(\max(0,U_i-l_i)\) |
| (b_i\Rightarrow r_i=q_i-d_i) | \(\max(0,-L_i,U_i-l_i)\) em cada sentido |
| (c_i\Rightarrow d_i\ge q_i) | \(\max(0,u_i-L_i)\) |
| (c_i\Rightarrow r_i=0) | \(u_i\) em cada sentido |

Esses valores cobrem os **outros** regimes quando a implicação está desligada: em todos os regimes (0\le r_i\le q_i\le u_i); no regime (a), (r_i-q_i+d_i=d_i); no regime (c), (r_i-q_i+d_i=d_i-q_i). `deficit_bounds` servem exclusivamente para calcular constantes; nenhuma cota é aplicada à variável (d_i). Nenhuma cota de `realized_bounds` é aplicada a (r_i), que é criada sem limites próprios. O clipping e o domínio público implicam todas as propriedades usadas acima.

## Objetivo e limite

\[
FOPT=\beta_0+\sum_{i=1}^{8}\beta_i real_i+\beta_{Z2}Z2,
\qquad FWIT=\gamma_0+\gamma_ZZ.
\]

Os coeficientes são lidos de `layer2`. Maximiza-se `FOPT` sujeito a `FWIT <= epsilon_FWIT`, passado por `--epsilon`. WBHP pertence apenas à validação física externa e não aparece no MILP. Não são usadas soluções de referência nem resultados OPM.
