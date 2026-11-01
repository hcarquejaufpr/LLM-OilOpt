# Formulação MILP — P2 público v3

Fonte única dos parâmetros: `shared/problema_p2_public.json`. O executável lê esse arquivo em cada execução. Usa `scipy.optimize.milp`, com o resolvedor HiGHS.

## Variáveis e restrições

Para cada poço, `q_i` tem exatamente o domínio público; `q_1=79.5`. Impõe-se `Σq_i≤636`. As igualdades das features são `20Q=Σq_i−600`, `20own_i=q_i−80`, `R=Σr_i`, `20Z=R−600` e `20real_i=r_i−80`.

O déficit `d_i` é a combinação linear dos coeficientes de `layer1_deficit_models`. Nos poços sem coeficientes, `d_i=0` e `r_i=q_i`. Em INJECT3, os termos quadráticos usam os valores públicos interpolados de `own2_INJECT3` e `Q2`.

Para cada função PWL de nós `(x_j,y_j)`, são criados pesos contínuos `λ_j≥0` e seletores binários `b_j` para os segmentos. As equações são `Σλ_j=1`, `Σb_j=1`, `x=Σx_jλ_j`, `y=Σy_jλ_j` e `λ_j≤b_(j−1)+b_j`, com seletores ausentes nas extremidades iguais a zero. Assim, apenas os dois nós adjacentes do segmento escolhido podem ter peso positivo. As três funções preservam literalmente nós e valores do JSON público.

## Clipping exato

Para os poços com déficit modelado, há seletores binários `z0+z1+z2=1`:

| Seletor | Condição | Vazão realizada |
| --- | --- | --- |
| `z0` | `d≤0` | `r=q` |
| `z1` | `0≤d≤q` | `r=q−d` |
| `z2` | `d≥q` | `r=0` |

Nos limites compartilhados, as fórmulas coincidem. Logo, `r=q−clip(d,0,q)` em todo ponto factível.

Sejam `[L,U]=deficit_bounds_i` e `[a,b]` o domínio público de `q_i`. Os Big-M de cada implicação são:

| Implicação | Desativação suficiente |
| --- | --- |
| `z0⇒d≤0` | `max(0,U)` |
| `z1⇒d≥0` | `−min(0,L)` |
| `z1⇒d−q≤0` | `max(0,U−a)` |
| `z2⇒d−q≥0` | `max(0,b−L)` |
| `z0⇒r−q=0` | `b` para ambas as faces |
| `z1⇒r−q+d=0` | `max(|L|,|U|,b+max(0,−L))` para ambas as faces |
| `z2⇒r=0` | `b` para ambas as faces |

As duas faces das igualdades são ativadas com Big-M simétrico. Esses valores cobrem os valores das expressões nas outras regiões ativas, mesmo sem limites explícitos em `d` ou `r`. `deficit_bounds` entram apenas no cálculo das constantes; `realized_bounds` não são necessários nesta escolha de Big-M. Nenhum dos dois é imposto como restrição adicional. As implicações ativas, os domínios de `q` e as igualdades públicas definem o conjunto factível.

## Objetivo

`FOPT=intercept+Σ beta_i·real_i+beta_Z2·Z2`, com os coeficientes de `layer2.FOPT`; `FWIT=intercept+beta_Z·Z`, com os coeficientes de `layer2.FWIT`. Maximiza-se `FOPT`, sujeito a `FWIT≤epsilon_FWIT` fornecido na linha de comando. WBHP é apenas uma validação física externa e não entra no MILP.

Após a solução, `solve.py` recalcula déficit, clipping, `FOPT` e `FWIT` diretamente das vazões e dos dados públicos. Divergências acima da tolerância numérica causam erro. A saída contém as variáveis de decisão, valores realizados, features relevantes, déficits, objetivo, epsilon e gap do resolvedor.
