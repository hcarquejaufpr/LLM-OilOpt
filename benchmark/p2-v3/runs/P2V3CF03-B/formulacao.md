# Formulação MILP — P2 público v3

Fonte única dos dados: `shared/problema_p2_public.json`. O programa lê esse arquivo em cada execução. As variáveis de vazão estão na ordem `wells` do arquivo público.

## Domínio e features

`q_INJECT1 = 79.5`; para os demais poços, `q_i` usa exatamente os intervalos em `decision_variables.domains`. Impõe-se `sum_i q_i <= 636`. As igualdades de normalização são `20 Q = sum_i q_i - 600` e `20 own_i = q_i - 80`.

Para cada poço, `deficit_i` é uma igualdade afim com os coeficientes de `layer1_deficit_models`. Os modelos vazios dão déficit zero. Em INJECT3 entram ainda `own2_INJECT3` e `Q2`, interpolados nos nós públicos.

## PWL exata nos nós públicos

Para cada função, sejam `x_j,v_j` os nós e valores públicos, `lambda_j >= 0` e `s_k` binárias dos segmentos. Impõe-se `sum_j lambda_j = 1`, `sum_k s_k = 1`, `x = sum_j x_j lambda_j`, `v = sum_j v_j lambda_j` e `lambda_j <= s_(j-1)+s_j` (com termos de borda omitidos). Logo, somente os dois nós de um segmento ativo podem ter peso positivo. Isso implementa a interpolação linear adjacente para `own2_INJECT3`, `Q2` e `Z2`.

## Clipping exato em três regimes

Para cada poço, usam-se binárias `n_i,m_i,h_i` com soma 1. Os regimes são:

| Binária ativa | Condição | Vazão realizada |
| --- | --- | --- |
| `n_i` | `d_i <= 0` | `r_i = q_i` |
| `m_i` | `0 <= d_i <= q_i` | `r_i = q_i - d_i` |
| `h_i` | `d_i >= q_i` | `r_i = 0` |

Os empates nas fronteiras podem ativar qualquer regime adjacente, com o mesmo valor de `r_i`. `0 <= r_i <= q_i` é consequência da definição de clipping e fortalece a relaxação.

Se `[L_d,U_d]` vem de `deficit_bounds` e `[L_q,U_q]` é o domínio público de `q_i`, as constantes suficientes, armazenadas em `model.json` e calculadas em `solve.py`, são:

| Implicação | Desativação Big-M |
| --- | --- |
| `n: d <= 0` | `max(0,U_d)` |
| `m: d >= 0` | `max(0,-L_d)` |
| `m: d <= q` | `max(0,U_d-L_q)` |
| `h: d >= q` | `max(0,U_q-L_d)` |
| `n: r = q` | `U_q` para `r-q >= -M`; `r <= q` já é global |
| `m: r = q-d` | `max(0,U_q-L_d)` para `r-q+d >= -M`; `max(0,U_d)` para `r-q+d <= M` |
| `h: r = 0` | `U_q` para `r <= M` |

Todas essas constantes são válidas por aritmética intervalar usando `0 <= r <= q`, os domínios de `q` e os limites auxiliares de déficit. `deficit_bounds` não são cotas da variável `deficit_i` nem restrições. `realized_bounds` tampouco são cotas ou restrições de `r_i`. O limite de variável `0 <= r_i <= U_q` é consequência direta do clipping; o limite superior mais forte vem de `r_i <= q_i`.

## Camada realizada e otimização

`R = sum_i r_i`, `20 Z = R - 600` e `20 real_i = r_i - 80`. A função `Z2` interpola os nós públicos no argumento `Z`. `FOPT = intercept + sum_i beta_ri real_i + beta_Z2 Z2`, com os coeficientes de `layer2.FOPT`. `FWIT = intercept + beta_Z Z`, com os coeficientes de `layer2.FWIT`. Maximiza-se `FOPT` sujeito a `FWIT <= epsilon_FWIT`, fornecido em `--epsilon`.

WBHP é apenas validação física externa e não aparece no MILP. Nenhum resultado de referência ou OPM é usado. O arquivo de saída calcula novamente déficits, clipping, FOPT e FWIT a partir das vazões da solução e falha se houver divergência material.
