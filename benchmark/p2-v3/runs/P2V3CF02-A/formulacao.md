# Formulação P2 v3 — injeção de água

Fonte: `shared/problema_p2_public.json`. A implementação usa MILP com segmentos adjacentes para cada função PWL e maximiza `FOPT` sob o orçamento `epsilon_FWIT` recebido na linha de comando. O horizonte de 3600 dias já está incorporado nos coeficientes públicos de FOPT e FWIT; não há multiplicação adicional.

## Variáveis e domínio

Para cada poço `i`, `q_i` é a vazão solicitada, `d_i` o déficit do modelo de primeira camada e `r_i` a vazão realizada. `q_1=79.5`; os demais `q_i` pertencem aos intervalos em `decision_variables.domains`. Impõe-se `sum_i q_i <= 636`. O modelo também usa as igualdades

```
Q     = (sum_i q_i - 600)/20
own_i = (q_i - 80)/20
Z     = (sum_i r_i - 600)/20
real_i= (r_i - 80)/20
```

Somente `own_3` aparece com coeficiente não nulo no déficit. Os déficits são exatamente as combinações de `layer1_deficit_models`, com campos ausentes iguais a zero:

```
d_i = intercept_i + beta_own_i*own_i + beta_Q_i*Q
    + beta_own2_i*PWL_own2(own_3) + beta_Q2_i*PWL_Q2(Q).
```

Em particular, `d_2=d_4=d_5=0`.

## Funções PWL

Para cada tabela pública com nós `(x_j,y_j)`, introduzimos pesos `lambda_j >= 0` e binárias de segmento `s_k`, `k=0,...,n-2`:

```
sum_j lambda_j = 1,  sum_k s_k = 1
x = sum_j x_j lambda_j,  y = sum_j y_j lambda_j
lambda_0 <= s_0, lambda_{n-1} <= s_{n-2}
lambda_j <= s_{j-1}+s_j,  j=1,...,n-2.
```

Isso permite somente os dois nós adjacentes de um segmento. As tabelas são `own2_INJECT3`, `Q2` e `Z2`; seus valores são interpolados linearmente entre os nós públicos. Os intervalos dos nós cobrem os argumentos alcançáveis pelo conjunto factível: `own_3` vai de `-1.025` a `3.0125`; `Q` vai de `-8.95` a `1.8`; o limite inferior conservador de `Z`, calculado por `sum_i(q_i^min-max(d_i^max,0))`, coincide com o primeiro nó de `Z2` à precisão dos dados, e `Z<=1.8` decorre de `r_i<=q_i` e da capacidade. A tabela PWL define a função no intervalo necessário; não se usam `realized_bounds` como limites de variáveis.

## Clipping exato

Para cada poço, três binárias `a_i,b_i,c_i` somam 1 e representam, respectivamente, `d_i<=0`, `0<=d_i<=q_i`, `d_i>=q_i`. A implicação ativa é:

```
a_i=1  =>  d_i<=0,            r_i=q_i
b_i=1  =>  0<=d_i<=q_i,      r_i=q_i-d_i
c_i=1  =>  d_i>=q_i,         r_i=0.
```

Nas fronteiras, mais de um regime pode representar o mesmo ponto, sem alterar `r_i=q_i-clip(d_i,0,q_i)`. Também se impõe `0<=r_i<=q_i`, consequência direta dessa identidade. O limite de variável `r_i<=q_i^max` é redundante e decorre do domínio de `q_i`; ele não é `realized_bounds`.

As implicações são linearizadas com as seguintes desigualdades. Se `L_i,U_i` são os **limites auxiliares** em `deficit_bounds` e `l_i,u_i` são os extremos do domínio de `q_i`, definem-se:

```
A=max(U_i,0), B=max(-L_i,0), C=max(U_i-l_i,0),
D=max(u_i-L_i,0), E=max(U_i,0).

d_i                    <= A(1-a_i)
r_i-q_i                >= -u_i(1-a_i)
-d_i                   <= B(1-b_i)
d_i-q_i                <= C(1-b_i)
r_i-q_i+d_i            >= -D(1-b_i)
r_i-q_i+d_i            <= E(1-b_i)
q_i-d_i                <= D(1-c_i)
r_i                    <= u_i(1-c_i).
```

Cada constante cobre o extremo da expressão correspondente no regime inativo, usando `L_i<=d_i<=U_i` **somente para certificar o Big M**, `l_i<=q_i<=u_i` e `0<=r_i<=q_i`. Nenhuma desigualdade `L_i<=d_i<=U_i` entra no MILP. Também não entram limites de `realized_bounds`. Os limites auxiliares são declarados válidos no arquivo público para o conjunto factível completo; a auditoria numérica confere as soluções obtidas, sem substituir essa garantia global.

## Objetivo e orçamento

Os coeficientes públicos da segunda camada definem:

```
FOPT = alpha_F + sum_{i=1}^8 beta_i*(r_i-80)/20
       + beta_Z2*PWL_Z2(Z)
FWIT = alpha_W + beta_Z*Z
maximize FOPT
FWIT <= epsilon_FWIT.
```

`WBHP` não é variável nem restrição: seu limite é destinado à validação física externa. Não há outras restrições de viabilidade além das igualdades definidoras e do conjunto factível público. Em particular, `deficit_bounds` e `realized_bounds` não reduzem o conjunto factível.

## Execução e verificação

`python3 solve.py --epsilon <valor> --output <arquivo.json>` usa `scipy.optimize.milp` com HiGHS. O JSON de saída contém status, vazões, soma realizada `R`, métricas, gap e resíduos de uma recomputação independente das equações PWL, dos déficits e do clipping. `evidence.json` resume as cinco execuções públicas; as soluções específicas não são codificadas no solver.
