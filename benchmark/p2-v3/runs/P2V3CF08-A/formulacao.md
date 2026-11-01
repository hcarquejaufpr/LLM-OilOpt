# Formulação MILP/PWL — P2 público v3

## Dados e objetivo

A instância é lida de `shared/problema_p2_public.json`. Para cada poço (i), (q_i) é a vazão solicitada, (d_i) é o déficit previsto e (r_i) é a vazão realizada. (q_1=79{,}5); os demais (q_i) pertencem aos intervalos públicos. A única restrição conjunta de vazão solicitada é (sum_i q_i\le636). Defina

\[
Q=(\sum_i q_i-600)/20,\quad o_i=(q_i-80)/20,\quad
R=\sum_i r_i,\quad Z=(R-600)/20,\quad t_i=(r_i-80)/20.
\]

Os coeficientes abaixo são lidos do JSON, sem valores de solução embutidos. Para cada poço,

\[
d_i=a_i+b_i o_i+c_iQ+e_i\operatorname{PWL}_{o_3^2}(o_3)
       +f_i\operatorname{PWL}_{Q^2}(Q),
\]

com coeficientes ausentes iguais a zero. As três funções PWL usam *exatamente* os nós e valores públicos e interpolação linear entre nós adjacentes. Os nomes `own2` e `Q2` identificam essas tabelas; o MILP não calcula quadrados analíticos entre os nós.

\[
r_i=q_i-\operatorname{clip}(d_i,0,q_i),
\]
\[
\mathrm{FOPT}=a_F+\sum_{i=1}^8\beta_i t_i
              +\beta_{Z2}\operatorname{PWL}_{Z^2}(Z),
\qquad
\mathrm{FWIT}=a_W+\beta_Z Z.
\]

O modelo maximiza FOPT sujeito a (\mathrm{FWIT}\le\epsilon_{\mathrm{FWIT}}). `epsilon_FWIT` é argumento de execução, sem restrição aos cinco valores da lista. WBHP não participa do MILP.

## Interpolação PWL com adjacência

Para uma tabela ((x_k,v_k)), usam-se (\lambda_k\ge0) e indicadores binários (s_j), um por segmento. Impõem-se

\[
\sum_k\lambda_k=1,\quad \sum_j s_j=1,\quad
x=\sum_k x_k\lambda_k,\quad y=\sum_k v_k\lambda_k,
\quad \lambda_k\le s_{k-1}+s_k,
\]

omitindo os segmentos inexistentes nas extremidades. Assim, somente os dois nós de um mesmo segmento podem ter peso positivo. Isso equivale a SOS2 e determina exatamente a interpolação pública. O domínio dos nós inclui os valores possíveis dos argumentos conforme os dados públicos.

## Clipping em três regimes

Para cada poço com modelo de déficit não vazio, sejam (z_i^-,z_i^0,z_i^+\in\{0,1\}), com soma 1. As implicações são

| Indicador | Condição | Vazão realizada |
|---|---|---|
| (z_i^-=1) | (d_i\le0) | (r_i=q_i) |
| (z_i^0=1) | (0\le d_i\le q_i) | (r_i=q_i-d_i) |
| (z_i^+=1) | (d_i\ge q_i) | (r_i=0) |

Na igualdade entre regimes, ambas as escolhas produzem o mesmo (r_i). Para os modelos de déficit vazios, (d_i=0) e a identidade é diretamente (r_i=q_i). Os limites de variável (0\le r_i\le\overline q_i) são consequências da identidade de clipping, não dos `realized_bounds` auxiliares.

Cada implicação (z=1\Rightarrow h\le0) é escrita (h\le M(1-z)); uma igualdade usa as duas desigualdades de sinais opostos. Para (q_i\in[a,b]) e os limites auxiliares válidos (d_i\in[L,U]), a implementação usa os seguintes limites para os respectivos lados inativos, acrescidos de uma margem numérica de (10^{-7}):

| Implicação ativa | Lado inativo limitado | (M) antes da margem |
|---|---|---:|
| inferior: (d\le0) | (d) | \(\max(0,U)\) |
| médio: (d\ge0) | \(-d\) | \(\max(0,-L)\) |
| médio: (d\le q) | \(d-q\) | \(\max(0,U-a)\) |
| superior: (d\ge q) | \(q-d\) | \(\max(0,b-L)\) |
| inferior: (r=q) | \(|r-q|\) | \(\min(b,\max(0,U))\) |
| médio: (r=q-d) | \(|r-q+d|\) | \(\max(0,-L,U-a)\) |
| superior: (r=0) | \(|r|\) | \(b\) |

Esses limites cobrem os outros regimes de clipping: (r-q=-\operatorname{clip}(d,0,q)); (r-q+d=d-\operatorname{clip}(d,0,q)), que é (d) no regime inferior, zero no médio e (d-q) no superior. A margem apenas amplia Big-M; quando o indicador vale 1, a implicação mantém a igualdade ou desigualdade exata. `deficit_bounds` e `realized_bounds` não são impostas como restrições de factibilidade. Os primeiros servem somente ao cálculo de Big-M; os segundos são apenas conferidos na auditoria da solução.

## Implementação e auditoria

`solve.py` monta uma matriz esparsa MILP e usa SciPy/HiGHS. `model.json` contém a matriz explícita da instância com o epsilon P1, os tipos de variáveis, os limites, a função objetivo e o hash SHA-256 do JSON público. Para outro epsilon, apenas o lado direito da linha `FWIT_limit` muda. As entradas `null` para limites na matriz significam infinito.

Execute `python3 solve.py --epsilon 1728051.611471 --output resultado.json`. A saída contém status, vazões, R, FOPT, FWIT, gap e resíduos calculados novamente a partir das equações públicas. `evidence.json` registra a execução dos cinco epsilons públicos. Pequenos resíduos são tolerâncias numéricas do solver, não restrições adicionais.
