# Contrato de entrega — P2 v3 A/B

As duas condições obedecem ao mesmo contrato.

## Arquivos obrigatórios

- `formulacao.md`
- `solve.py`
- `model.json`
- `evidence.json`

Uma entrega completa deve existir na raiz de trabalho ou em uma única subpasta inequívoca.

## Interface de solve.py

`python3 solve.py --epsilon <valor> --output <arquivo.json>`

O processo deve terminar com código 0 quando o problema é resolvido normalmente.

O arquivo de saída deve conter:
- `status`
- `epsilon_FWIT`
- `requested`
- `realized`
- `R`
- `FOPT`
- `FWIT`
- `solver_gap`

Campos adicionais são permitidos para auditoria.

## Regras

- não modificar `shared/problema_p2_public.json`;
- não buscar arquivos externos;
- não usar rede;
- não incorporar WBHP ao MILP;
- não impor `deficit_bounds` como restrição de factibilidade sobre `deficit_i` (uso permitido: derivar/verificar Big-M); o mesmo vale para `realized_bounds` sobre `r_i`;
- não codificar soluções P1–P5 conhecidas;
- usar a definição pública v3 de FOPT, inclusive `real_i=(r_i-80)/20`;
- a mesma formulação deve funcionar para todos os epsilon públicos;
- erros de formulação ou execução contam como resultado experimental e não justificam rerun.

## Arquivos copiados ao workspace de sessão (lista fechada)

- Ambas as condições: `shared/problema_p2_public.json` e o prompt/SDD da condição correspondente (`PROMPT_CONDICAO_A.md` ou `SDD_CONDICAO_B.md`), mais `CONTRATO_ENTREGA.md`.
- Somente condição B: `TDD_CONDICAO_B.md` e `B/public_tdd.py`.
- Nunca copiar: `evaluator_private/`, `relatorios_agentes/`, `GERAR_PROBLEMA_PUBLICO_V3.py`, `FORMULACAO_MATEMATICA_V3.md`, `HISTORICO_V1_V2_V3.md`, qualquer arquivo v1/v2, `auth/`.
