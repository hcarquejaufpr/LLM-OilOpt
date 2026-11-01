# Formulação MILP/PWL — P2 v3

## Dados e conjunto factível

Todos os números vêm de `shared/problema_p2_public.json`. Seja \(I=\{1,\ldots,8\}\), com a ordem dos poços `INJECT1` a `INJECT8`. As vazões solicitadas \(q_i\) obedecem aos domínios do arquivo: \(q_1=79{,}5\), limites individuais para \(i=2,\ldots,8\), e \(\sum_i q_i\le636\). O parâmetro \(\epsilon_{FWIT}\) é dado na chamada do programa. Fora das definições funcionais abaixo, as únicas restrições são esses domínios, a capacidade solicitada e \(FWIT\le\epsilon_{FWIT}\).

Os limites `deficit_bounds` e `realized_bounds` são auxiliares. Nenhum deles é imposto como restrição sobre \(d_i\) ou \(r_i\). O limite de WBHP é usado somente em validação física externa e não aparece no MILP.

## Features e déficits

Defina

\[
Q=\frac{\sum_iq_i-600}{20},\qquad
u_i=\frac{q_i-80}{20},\qquad
Z=\frac{\sum_ir_i-600}{20},\qquad
v_i=\frac{r_i-80}{20}.
\]

As variáveis auxiliares \(u_3=\operatorname{PWL}_{\mathrm{own2\_INJECT3}}(\nu_3)\), \(s_Q=\operatorname{PWL}_{Q2}(Q)\) e \(s_Z=\operatorname{PWL}_{Z2}(Z)\) usam exatamente os nós e valores públicos. Os déficits são

\[
\begin{aligned}
d_1 &= 1.1713666838726333+0.595529655325077Q,\\
d_2 &=0,\\
d_3 &=-0.015303143431095477+0.19911165548663023\nu_3
       -0.006664634656562543Q+0.17884392421295806u_3
       +0.019654062343378798s_Q,\\
d_4&=d_5=0,\\
d_6&=1.237674360935837+2.0895406180704765\nu_6+0.1880912450856571Q,\\
d_7&=0.5982398897849602+1.3233182936757202\nu_7-0.012674986258828236Q,\\
d_8&=1.6365876680135958+2.450483665840999\nu_8+0.15031172964269102Q.
\end{aligned}
\]

### PWL com segmentos adjacentes

Para cada função pública \(f\), com nós \(x_j\) e valores \(y_j\), o modelo cria pesos \(\lambda_j\in[0,1]\) e binárias \(b_k\) para os segmentos \([x_k,x_{k+1}]\):

\[
\sum_j\lambda_j=1,\quad \sum_k b_k=1,\quad
x=\sum_j x_j\lambda_j,\quad f(x)=\sum_j y_j\lambda_j,
\]
\[
\lambda_0\le b_0,\quad \lambda_{n-1}\le b_{n-2},\quad
\lambda_j\le b_{j-1}+b_j\quad(1\le j<n-1).
\]

Assim, só os dois extremos de um mesmo segmento podem ter peso positivo. Esta é uma formulação MILP equivalente à adjacência SOS2; a interpolação é exata para as funções PWL especificadas, sem aproximação quadrática adicional.

## Clipping exato em três regimes

Para cada poço, três binárias \(a_{i0},a_{i1},a_{i2}\) somam 1. Cada regime ativa as implicações:

| Regime | Condição do déficit | Vazão realizada |
|---|---|---|
| 0 | \(d_i\le0\) | \(r_i=q_i\) |
| 1 | \(0\le d_i\le q_i\) | \(r_i=q_i-d_i\) |
| 2 | \(d_i\ge q_i\) | \(r_i=0\) |

Nas fronteiras, dois regimes podem representar a mesma solução física. Em todos os casos, \(r_i=q_i-\operatorname{clip}(d_i,0,q_i)\).

Uma implicação ativa \(g\le0\) é escrita \(g\le M(1-a)\); uma igualdade \(g=0\) usa \(-M(1-a)\le g\le M(1-a)\). Para cada poço, sejam \([L_i,U_i]\) os `deficit_bounds` auxiliares e \([\ell_i,h_i]\) o domínio de \(q_i\). Os Big-M usados, na ordem das implicações acima, são:

| Implicação ativa | Big-M |
|---|---:|
| \(d_i\le0\) | \(\max(0,U_i)\) |
| \(r_i=q_i\) | \(h_i\) |
| \(d_i\ge0\) | \(\max(0,-L_i)\) |
| \(d_i\le q_i\) | \(\max(0,U_i-\ell_i)\) |
| \(r_i=q_i-d_i\) | \(\max(h_i-L_i,\;h_i-\ell_i+U_i)\) |
| \(d_i\ge q_i\) | \(\max(0,h_i-L_i)\) |
| \(r_i=0\) | \(h_i\) |

Essas constantes cobrem qualquer ramo inativo: em toda solução do clipping vale \(0\le r_i\le q_i\le h_i\), e os intervalos auxiliares cobrem os valores dos modelos de déficit. Essa observação serve apenas à prova dos Big-M, sem adicionar limites a \(r_i\) ou \(d_i\). Os valores numéricos por poço estão em `model.json`, campo `clipping_big_m`.

## Objetivo e orçamento

\[
FOPT=501555.5938340376+\sum_{i=1}^{8}\beta_i v_i
      -167.73716538798362s_Z,
\]

com \((\beta_1,\ldots,\beta_8)=(-14719.837617462217,235.1707114443161,1599.6347095184774,318.9085949158334,570.0842702497193,1606.473270241966,1703.0746679149643,1783.656679566027)\). Além disso,

\[
FWIT=2159936.0100392425+71980.73309472566Z,
\qquad FWIT\le\epsilon_{FWIT}.
\]

O objetivo é **maximizar \(FOPT\)**. O valor `R` no JSON de saída é \(\sum_i r_i\).

## Implementação e auditoria

`solve.py` lê o arquivo público a cada execução, monta a matriz esparsa e resolve com SciPy `milp` (HiGHS). Exemplo: `python3 solve.py --epsilon 1943993.810755 --output resultado.json`. Requer `scipy` e `numpy`. O `solver_gap` é o gap relativo informado pelo HiGHS. O arquivo `model.json` contém a matriz esparsa gerada para o epsilon P3; a linha `epsilon_FWIT` recebe o valor informado em cada nova execução. `evidence.json` contém as cinco execuções públicas e resíduos recalculados independentemente a partir das vazões e da interpolação dos nós.
