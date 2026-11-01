# Formulação MILP — P2 público v3

Fonte única dos dados: `shared/problema_p2_public.json`. O script carrega esse arquivo em tempo de execução. Todas as vazões são contínuas.

## Variáveis e conjunto factível

Para cada `INJECTi`, `q_i` tem exatamente o domínio em `decision_variables`; `q_1=79.5`. Impõe-se `sum_i q_i <= 636`. Definem-se

\[
Q=(\sum_i q_i-600)/20,\quad own_i=(q_i-80)/20,\quad
R=\sum_i r_i,\quad Z=(R-600)/20,\quad real_i=(r_i-80)/20.
\]

Cada déficit `d_i` é a expressão afim de `layer1_deficit_models`, usando `own2_INJECT3` e `Q2` somente onde seus coeficientes aparecem. Um modelo vazio define `d_i=0`.

## Funções lineares por partes

Para cada função pública (`own2_INJECT3`, `Q2`, `Z2`) e seus pares de nós `(x_j,y_j)`, usam-se pesos contínuos `lambda_j >= 0` e indicadores binários de segmento `s_k`. As equações são `sum_j lambda_j=1`, `sum_k s_k=1`, `x=sum_j x_j lambda_j`, `y=sum_j y_j lambda_j` e `lambda_j <= s_{j-1}+s_j`, omitindo índices inexistentes. Logo só dois nós adjacentes podem ter peso positivo; a interpolação é exata entre os nós fornecidos. Os arrays são lidos sem alterar valores de `pwl_functions`.

## Clipping exato em três regimes

Para cada poço, `t_i=clip(d_i,0,q_i)` e `r_i=q_i-t_i`. `t_i>=0` e `t_i<=q_i` são propriedades intrínsecas do clip. Indicadores binários `a_i,b_i,c_i` somam 1:

| Indicador ativo | Condições | Valor de `t_i` |
|---|---|---|
| `a_i` | `d_i<=0` | `0` |
| `b_i` | `0<=d_i<=q_i` | `d_i` |
| `c_i` | `d_i>=q_i` | `q_i` |

As constantes Big-M são calculadas para cada poço a partir de `L_i,U_i=deficit_bounds[i]` e `l_i,u_i=domain(q_i)`. Na ordem das oito implicações acima implementadas em `add_clip`, elas são:

| Implicação | Constante suficiente |
|---|---|
| `a => d<=0` | `max(U,0)` |
| `a => t=0` | `u` |
| `b => d>=0` | `max(-L,0)` |
| `b => d<=q` | `max(U-l,0)` |
| `b => t-d<=0` | `max(u-L,0)` |
| `b => d-t<=0` | `max(U,0)` |
| `c => d>=q` | `max(u-L,0)` |
| `c => t=q` | `u` |

Essas constantes tornam cada implicação inativa válida para qualquer valor nos intervalos auxiliares. `deficit_bounds` **não** são limites da variável `d_i`; `realized_bounds` **não** são limites da variável `r_i`. O código não cria restrições adicionais com esses dados. As únicas restrições de factibilidade são os domínios de `q`, a capacidade solicitada, as igualdades de definição, o clipping e o teto de FWIT. WBHP não entra no MILP.

## Objetivo e execução

`FOPT = intercept + sum_i beta_ri * real_i + beta_Z2 * Z2`, com todos os coeficientes de `layer2.FOPT`. `FWIT = intercept + beta_Z * Z`, com `layer2.FWIT`. Maximiza-se FOPT sujeito a `FWIT <= epsilon_FWIT`.

Execute `python3 solve.py --epsilon 1943993.810755 --output resultado.json`. O programa usa SciPy/HiGHS, exige prova de otimalidade para terminar com código 0 e grava as vazões, déficits, features, FOPT, FWIT, gap, Big-M e resíduo do clipping em JSON. Nenhuma solução conhecida ou resultado OPM é usado.
