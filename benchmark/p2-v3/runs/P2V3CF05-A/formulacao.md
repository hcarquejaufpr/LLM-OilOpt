# Formulação MILP/PWL — P2 público v3

## Dados e variáveis

Todos os coeficientes, domínios e nós PWL são lidos de `shared/problema_p2_public.json`. Para cada poço (i), (q_i) é a vazão solicitada, (r_i) a realizada e (d_i) o déficit calculado pela primeira camada. (q_1=79{,}5); os demais (q_i) pertencem aos respectivos intervalos públicos. Define-se

\[
Q=(\sum_i q_i-600)/20,\quad Z=(\sum_i r_i-600)/20,\quad
own_i=(q_i-80)/20,\quad real_i=(r_i-80)/20.
\]

O orçamento é \(FWIT\leq\epsilon_{FWIT}\), e a capacidade solicitada é \(\sum_iq_i\leq636\). O objetivo é maximizar FOPT.

## Primeira camada e PWL

Os déficits são exatamente as combinações lineares dos coeficientes públicos. Para poços com modelo vazio, \(d_i=0\). Os termos `own2_INJECT3`, `Q2` e `Z2` são as **interpolações lineares dos valores tabulados**, sem substituir seus valores por quadrados analíticos.

Para cada função PWL com nós \((x_j,y_j)\), usam-se \(\lambda_j\geq0\) e indicadores binários \(s_k\) de segmento:

\[
\sum_j\lambda_j=1,\quad \sum_ks_k=1,\quad
x=\sum_jx_j\lambda_j,\quad y=\sum_jy_j\lambda_j,\quad
\lambda_j\leq s_{j-1}+s_j,
\]

omitindo índices de segmento fora da tabela. Assim, somente os dois nós de um único segmento adjacente podem receber peso positivo. Os três pares \((x,y)\) são \((own_3,own2_3)\), \((Q,Q2)\) e \((Z,Z2)\).

## Clipping exato

Em cada poço, três indicadores \(a_i,b_i,c_i\in\{0,1\}\) satisfazem \(a_i+b_i+c_i=1\). As implicações são:

\[
\begin{array}{c|c|c}
\text{regime}&\text{condição}&\text{vazão realizada}\\\hline
a_i=1&d_i\leq0&r_i=q_i\\
b_i=1&0\leq d_i\leq q_i&r_i=q_i-d_i\\
c_i=1&d_i\geq q_i&r_i=0
\end{array}
\]

Nos limiares, dois regimes podem descrever o mesmo ponto, sem alterar \(r_i=q_i-\operatorname{clip}(d_i,0,q_i)\). O código lineariza cada lado de cada implicação com Big-M. Escrevendo \(q_L,q_U\) para o domínio de \(q_i\) e \(d_L,d_U\) para `deficit_bounds`, os maiores módulos usados quando uma implicação está inativa são:

| Implicação ativa | Expressão | Limite inativo usado |
|---|---|---|
| \(a=1\) | \(d\leq0\) | \(\max(0,d_U)\) |
| \(a=1\) | \(r-q\leq0\) | \(\max(0,q_U-q_L)\) |
| \(a=1\) | \(r-q\geq0\) | \(q_U\) |
| \(b=1\) | \(-d\leq0\) | \(\max(0,-d_L)\) |
| \(b=1\) | \(d-q\leq0\) | \(\max(0,d_U-q_L)\) |
| \(b=1\) | \(r-q+d\leq0\) | \(\max(0,q_U-q_L+d_U)\) |
| \(b=1\) | \(r-q+d\geq0\) | \(\max(0,q_U-d_L)\) |
| \(c=1\) | \(q-d\leq0\) | \(\max(0,q_U-d_L)\) |
| \(c=1\) | \(r\leq0\) | \(q_U\) |

Por exemplo, a primeira linha vira \(d+M a\leq M\), e uma igualdade \(e=0\) vira \(e+M_U z\leq M_U\) e \(e-M_L z\geq-M_L\). Os limites \(0\leq r_i\leq q_U\) usados para esses M são consequências da própria definição de clipping e de \(q_i\geq0\). Os intervalos `deficit_bounds` **não são restrições em** \(d_i\); `realized_bounds` **não são restrições em** \(r_i\). Nenhum deles reduz o conjunto factível. Os limites auxiliares do arquivo são assumidos válidos conforme sua definição pública; o código verifica sua ordem e usa seus extremos somente para construir os M.

## Segunda camada

Com os coeficientes públicos, a modelagem usa

\[
FOPT=\beta_0+\sum_{i=1}^8\beta_i\,real_i+\beta_{Z2}\,PWL_{Z2}(Z),
\qquad FWIT=\gamma_0+\gamma_Z Z.
\]

WBHP não pertence ao MILP; sua avaliação é validação física externa. Não se adiciona outra condição ao conjunto factível público.

## Execução

`python3 solve.py --epsilon 1943993.810755 --output resultado.json`

Requer Python 3, NumPy e SciPy com `scipy.optimize.milp` (HiGHS). O JSON de saída informa solução, lacuna do solver e variáveis intermediárias para auditoria. O arquivo `evidence.json` registra a execução dos cinco valores públicos e verificações independentes das equações.
