# P2 v3 — formulação MILP de injeção de água

## Dados e variáveis

Todos os coeficientes, domínios e nós PWL vêm de `shared/problema_p2_public.json`, reproduzido em `model.json`. Os poços são indexados por (i=1,\ldots,8). (q_i) é a vazão solicitada, (r_i) a realizada e (d_i) o déficit calculado pela primeira camada. (R=\sum_i r_i). As variáveis (Q), (Z), `own_i` e `real_i` representam exatamente

\[
Q=(\sum_i q_i-600)/20,\quad Z=(R-600)/20,\quad
\operatorname{own}_i=(q_i-80)/20,\quad
\operatorname{real}_i=(r_i-80)/20.
\]

Impõem-se (q_1=79.5), os sete domínios públicos de (q_2,\ldots,q_8) e (\sum_i q_i\le636). Não há outras restrições de operação.

## Primeira camada e PWL

Os déficits são as seguintes igualdades, usando os coeficientes públicos sem arredondamento no código:

\[
\begin{aligned}
d_1&=1.1713666838726333+0.595529655325077Q,\\
d_2&=0,\\
d_3&=-0.015303143431095477+0.19911165548663023\operatorname{own}_3
-0.006664634656562543Q+0.17884392421295806\operatorname{PWL}_{\mathrm{own2}}(\operatorname{own}_3)
+0.019654062343378798\operatorname{PWL}_{Q2}(Q),\\
d_4&=d_5=0,\\
d_6&=1.237674360935837+2.0895406180704765\operatorname{own}_6+0.1880912450856571Q,\\
d_7&=0.5982398897849602+1.3233182936757202\operatorname{own}_7-0.012674986258828236Q,\\
d_8&=1.6365876680135958+2.450483665840999\operatorname{own}_8+0.15031172964269102Q.
\end{aligned}
\]

Cada PWL usa exatamente os nós e valores publicados, com interpolação linear. Para uma PWL de (n+1) nós ((x_j,y_j)), criam-se pesos (\lambda_j\ge0) e binárias de segmento (s_k), (k=0,\ldots,n-1):

\[
\sum_j\lambda_j=1,\quad x=\sum_jx_j\lambda_j,\quad
y=\sum_jy_j\lambda_j,\quad \sum_ks_k=1,\quad
\lambda_j\le s_{j-1}+s_j,
\]

omitindo os segmentos inexistentes nas extremidades. Isso permite peso positivo somente em dois nós adjacentes do segmento escolhido, uma formulação MILP equivalente a SOS2. As PWL usadas são `own2_INJECT3`, `Q2` e `Z2`.

## Clipping exato

Para cada poço, sejam (q_i\in[\ell_i,u_i]) e os limites auxiliares válidos de déficit (L_i\le d_i\le U_i). Estes últimos **não são restrições** do modelo: servem somente para calcular Big-M. Três binárias (b_i^-,b_i^0,b_i^+\) somam 1 e selecionam, respectivamente,

\[
\begin{array}{c|c|c}
\text{regime}&\text{condição}&\text{vazão realizada}\\\hline
-&d_i\le0&r_i=q_i\\
0&0\le d_i\le q_i&r_i=q_i-d_i\\
+&d_i\ge q_i&r_i=0
\end{array}
\]

As igualdades dos regimes coincidem nas fronteiras. A implementação usa (0\le r_i\le u_i) e (r_i\le q_i), que são consequências do próprio clipping. As desigualdades condicionais implementadas para cada (i), com (b^-+b^0+b^+=1), são:

\[
\begin{aligned}
d&\le \max(U,0)(1-b^-),&q-r&\le u(1-b^-),\\
d&\ge \min(L,0)(1-b^0),&d-q&\le\max(U-\ell,0)(1-b^0),\\
r-q+d&\le\max(U,0)(1-b^0),&r-q+d&\ge-\max(u-L,0)(1-b^0),\\
d-q&\ge-\max(u-L,0)(1-b^+),&r&\le u(1-b^+).
\end{aligned}
\]

Cada Big-M é válido em todos os regimes desativados pelas faixas auxiliares públicas e pelo domínio de (q_i). Os valores numéricos estão em `model.json`. Os `deficit_bounds` não são impostos sobre (d_i). Os `realized_bounds` não são impostos sobre (r_i) nem usados como novos limites de viabilidade.

## Segunda camada e otimização

Com os coeficientes públicos (a_i) e (a_{Z2}):

\[
\begin{aligned}
\mathrm{FOPT}&=501555.5938340376+\sum_{i=1}^{8}a_i\operatorname{real}_i
-167.73716538798362\operatorname{PWL}_{Z2}(Z),\\
\mathrm{FWIT}&=2159936.0100392425+71980.73309472566Z.
\end{aligned}
\]

O vetor (a_i), na ordem dos poços, é

\[
(-14719.837617462217,\ 235.1707114443161,\ 1599.6347095184774,\ 318.9085949158334,\ 570.0842702497193,\ 1606.473270241966,\ 1703.0746679149643,\ 1783.656679566027).
\]

Maximiza-se FOPT sujeita a (\mathrm{FWIT}\le\epsilon_{\mathrm{FWIT}}). `epsilon_FWIT` é passado por `--epsilon`. WBHP é apenas uma validação física externa e não integra o MILP.

## Execução e auditoria

Requer Python 3 com NumPy e SciPy. Execute `python3 solve.py --epsilon 1943993.810755 --output resultado.json`. O programa lê o JSON público, monta o MILP e grava status, vazões, (R), FOPT, FWIT e gap do HiGHS. `evidence.json` registra as cinco execuções e resíduos de uma reavaliação independente das equações públicas a partir das vazões solicitadas. O status `optimal` é o certificado de otimalidade do solver dentro de `mip_rel_gap=1e-8` e das tolerâncias numéricas do HiGHS.
