# Formulação MILP do problema P2 v3

## Dados e conjunto factível

Todos os coeficientes, nós e valores PWL são lidos de `shared/problema_p2_public.json`. Para cada poço (i), (q_i) é a vazão solicitada, (d_i) o déficit calculado e (r_i) a vazão realizada. (q_1=79{,}5); os demais (q_i) obedecem aos intervalos públicos. Também vale (sum_i q_i\le636).

As features são

\[
Q=(\sum_i q_i-600)/20,\quad own_i=(q_i-80)/20,
\quad Z=(\sum_i r_i-600)/20,\quad real_i=(r_i-80)/20.
\]

O déficit de cada poço é exatamente a soma do intercepto e dos termos publicados para esse poço. Coeficientes ausentes valem zero. Assim, por exemplo,

\[
d_3=a_3+b_{3,o}own_3+b_{3,Q}Q
       +b_{3,o2}\operatorname{PWL}_{own2}(own_3)
       +b_{3,Q2}\operatorname{PWL}_{Q2}(Q).
\]

Os outros poços usam os termos `intercept`, `own` e `Q` presentes em seus respectivos objetos. Para `INJECT2`, `INJECT4` e `INJECT5`, isso dá (d_i=0).

## PWL adjacente

Para cada função PWL com nós (x_0<\cdots<x_n) e valores (v_0,\ldots,v_n), introduzimos, para cada segmento (s=0,\ldots,n-1), um binário (y_s) e pesos (a_s,b_s\ge0):

\[
\sum_s y_s=1,\quad a_s+b_s=y_s,\quad
x=\sum_s(a_s x_s+b_s x_{s+1}),\quad
p=\sum_s(a_s v_s+b_s v_{s+1}).
\]

Só os dois extremos de um segmento podem ter peso positivo. Isso equivale à adjacência SOS2, inclusive nos nós. A construção é aplicada a `own2_INJECT3`, `Q2` e `Z2`.

Os argumentos ficam cobertos pelos nós públicos: (own_3\in[-1{,}025,3{,}0125]), (Q\in[-8{,}95,1{,}8]). Para (Z), o limite superior decorre de (r_i\le q_i) e da capacidade. Um limite inferior conservador pode ser provado sem restringir (Z): se (Q\le-8), então (d_1<0) e (sum_i r_i\ge421-(U_3+U_6+U_7+U_8)>400{,}5); se (Q>-8), então (sum_i q_i>440) e (sum_i r_i>440-\sum_i\max(U_i,0)>417). Portanto (Z>-10), acima do primeiro nó de `Z2`.

## Clipping em três regimes

Defina (L_i=\texttt{deficit_bounds}[i][0]-10^{-7}), (U_i=\texttt{deficit_bounds}[i][1]+10^{-7}), e (\ell_i,u_i) como extremos do domínio público de (q_i). (L_i,U_i) só dimensionam constantes Big-M; **não são limites impostos à variável (d_i)**. Os limites de `realized_bounds` também **não são impostos a (r_i)**. Usamos apenas (0\le r_i\le q_i), que decorre da própria identidade de clipping.

Os binários (z_i^-,z_i^0,z_i^+\) somam 1 e representam, respectivamente:

| Regime | Condições ativas | Vazão realizada |
|---|---|---|
| baixo | (d_i\le0) | (r_i=q_i) |
| intermediário | (0\le d_i\le q_i) | (r_i=q_i-d_i) |
| alto | (d_i\ge q_i) | (r_i=0) |

As desigualdades exatas implementadas são, omitindo o índice (i):

\[
\begin{aligned}
d&\le U(1-z^-), & r&\ge q-u(1-z^-),\\
d&\ge L(1-z^0), & d-q&\le\max(0,U-\ell)(1-z^0),\\
r-q+d&\ge-(u-L)(1-z^0), & r-q+d&\le U(1-z^0),\\
d-q&\ge-(u-L)(1-z^+), & r&\le u(1-z^+).
\end{aligned}
\]

Junto a (0\le r\le q), cada regime ativo impõe exatamente sua linha da tabela. Os Big-M são válidos para os demais regimes pelas faixas auxiliares de déficit e pelo domínio de (q), sem mudar o conjunto factível. O termo (10^{-7}) protege contra arredondamento decimal dos limites auxiliares.

## Produção, orçamento e objetivo

Com os coeficientes públicos (c_i,c_{Z2},\alpha_i,\alpha_Z),

\[
FOPT=\alpha_0+\sum_{i=1}^8\alpha_i\frac{r_i-80}{20}
       +\alpha_{Z2}\operatorname{PWL}_{Z2}(Z),\qquad
FWIT=c_0+c_Z Z.
\]

Maximizamos `FOPT` sob `FWIT <= epsilon_FWIT`. A implementação minimiza `-FOPT`, de modo equivalente. `epsilon_FWIT` é um parâmetro da linha de comando; não há soluções ou políticas pré-codificadas. WBHP é apenas validação física externa e não entra no MILP.

## Auditoria numérica

`solve.py` usa `scipy.optimize.milp` com HiGHS. Após a solução, calcula independentemente déficits, clipping, PWL, FOPT e FWIT a partir das vazões solicitadas, registra os dois valores de objetivo/orçamento (avaliação direta e MILP) e o maior resíduo de vazão realizada. Um resíduo de limite de até (2\times10^{-7}) devolvido pelo solver é ajustado ao respectivo extremo do domínio antes da avaliação direta.
