# Rubrica semântica cega — P2 v3

Aplicação pós-campanha, sem abrir o mapa candidato→condição, avaliações privadas, timings ou tokens.

Uma entrega recebe `semantic_pass=PASS` somente se não houver falha essencial em:

1. arquitetura `requested -> realized -> FOPT/FWIT`;
2. domínios, `INJECT1=79.5`, capacidade e `FWIT<=epsilon`;
3. normalizações públicas, incluindo `real_i=(r_i-80)/20`;
4. déficits, clipping exato, PWL/SOS2 ou equivalente;
5. FOPT, FWIT e objetivo;
6. ausência de WBHP, hard-code ou informação privada;
7. `deficit_bounds` e `realized_bounds` somente como auxiliares, nunca como constraints.

Cada ficha usa os critérios C01–C23 da rubrica P2 v2 e acrescenta:

- C24: `deficit_bounds` não são constraints;
- C25: `realized_bounds` não são constraints.

`execution_pass` é verificado somente contra a interface e invariantes públicas, sem referência privada. `numeric_reference_pass` permanece `NOT_AVAILABLE` até autorização explícita de unblinding/M1–M4.
