# Formulação MILP — P2 público v3

Fonte única dos dados: `shared/problema_p2_public.json`. O programa lê esse arquivo em cada execução. `model.json` registra uma cópia dos dados, a estrutura do modelo e os Big-Ms calculados para auditoria.

## Variáveis e equações

Para cada poço (i), (q_i) tem exatamente o domínio público, com (q_1=79{,}5). Impomos (sum_i q_i\le636). As variáveis contínuas auxiliares obedecem às igualdades

\[
Q=(\sum_iq_i-600)/20,\quad own_i=(q_i-80)/20,\quad
R=\sum_ir_i,\quad Z=(R-600)/20,\quad real_i=(r_i-80)/20.
\]

Cada déficit (d_i) é igual à combinação afim de `intercept`, `own` e `Q` de `layer1_deficit_models`. Para INJECT3, acrescentamos exatamente os termos públicos `own2*PWL(own_3)` e `Q2*PWL(Q)`; modelos vazios dão (d_i=0).

Para cada função linear por partes, (lambda_j\ge0), (sum_j\lambda_j=1), (x=\sum_j x_j\lambda_j) e (y=\sum_j y_j\lambda_j). Binárias (s_k) selecionam um segmento: (sum_k s_k=1) e (lambda_j\le s_{j-1}+s_j), omitindo índices inexistentes nas extremidades. Assim só dois nós adjacentes podem ter peso positivo. Os nós e valores são lidos diretamente das três tabelas públicas.

## Clipagem em três regimes

Para cada poço, as binárias (b_i^L,b_i^M,b_i^H) somam 1. As implicações são:

| Regime | Condições quando ativo |
| --- | --- |
| (L) | (d_i\le0,\quad r_i=q_i) |
| (M) | (0\le d_i\le q_i,\quad r_i=q_i-d_i) |
| (H) | (d_i\ge q_i,\quad r_i=0) |

Nas fronteiras, regimes vizinhos podem coincidir, mas ambos produzem o mesmo (r_i=q_i-\operatorname{clip}(d_i,0,q_i)). Cada desigualdade (a\le0) é implementada como (a\le M(1-b)). Cada igualdade (a=0) usa (-M(1-b)\le a\le M(1-b)).

Se (d_i\in[d_L,d_U]) e (q_i\in[q_L,q_U]), os sete valores suficientes usados por `solve.py` são:

| Implicação | Big-M |
| --- | --- |
| (d\le0) | \(\max(0,d_U)\) |
| (d\ge0) | \(\max(0,-d_L)\) |
| (d\le q) | \(\max(0,d_U-q_L)\) |
| (d\ge q) | \(\max(0,q_U-d_L)\) |
| (r=q) | \(\max(0,d_U)\) |
| (r=q-d) | \(\max(0,-d_L,d_U-q_L)\) |
| (r=0) | \(q_U\) |

Para as três igualdades, os limites consideram o valor de (r) produzido pelos outros regimes: (q-r=\operatorname{clip}(d,0,q)), (r-q+d) vale (d) no regime baixo e (d-q) no alto, e (0\le r\le q_U). `interval_bounds` confere por aritmética intervalar que os limites públicos de déficit cobrem cada expressão em todo o domínio de (q), considerando capacidade e extremos das funções lineares por partes. `deficit_bounds` e `realized_bounds` não viram limites de variáveis nem linhas de restrição. O código não usa `realized_bounds`; não há restrição de WBHP.

## Camada 2 e objetivo

As igualdades para FOPT e FWIT usam os coeficientes de `layer2`: (FOPT=\beta_0+\sum_i\beta_i real_i+\beta_{Z2}PWL(Z)) e (FWIT=\gamma_0+\gamma_ZZ). O objetivo é maximizar FOPT sob (FWIT\le\epsilon_{FWIT}), informado por `--epsilon`. O arquivo de saída recalcula déficits, vazões realizadas, FOPT e FWIT de modo independente a partir da solução (q), e registra diferenças em relação às variáveis do MILP, status e gap do HiGHS.

Execução: `python3 solve.py --epsilon 1943993.810755 --output resultado.json`.
