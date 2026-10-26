# Seleção de artefatos para o pacote do artigo

**Data da triagem:** 2026-10-26  
**Origem:** `hcarquejaufpr/Pesquisa_Mestrado`, branch `main` (privado).  
**Situação:** seleção documental preliminar; nenhum arquivo experimental original foi copiado.

## Decisão de escopo

Separar o benchmark **P1** (MILP-PWL aditiva, quatro poços livres) da ampliação **P2** (sete variáveis livres, objetivos/restrições óleo–água). A campanha confirmatória do P1 tem protocolo e resultados; o P2 possui materiais de desenvolvimento que não devem ser misturados a P1.

A expressão “8×8” deve ser conferida antes da publicação: o problema de referência P1 usa oito poços injetores, mas somente quatro controles livres. Não rotular o benchmark P1 como oito variáveis livres.

## Arquivos candidatos (caminhos relativos ao repositório de origem)

| Prioridade | Origem | Destino proposto | Finalidade | Estado |
|---|---|---|---|---|
| 1 | `05_Experimentos_LLM/CAMPANHA_CONFIRMATORIA_20261019/runs/CF01-A/shared/problema.json` | `benchmark/p1/problem.json` | Instância matemática P1 | Selecionado; verificar proveniência e hash antes de copiar |
| 1 | `05_Experimentos_LLM/CAMPANHA_CONFIRMATORIA_20261019/PROTOCOLO_CONFIRMATORIO.md` | `docs/p1_confirmatory_protocol.md` | Protocolo A×B congelado | Selecionado; revisão editorial |
| 1 | `05_Experimentos_LLM/CAMPANHA_CONFIRMATORIA_20261019/CONGELAMENTO_CONFIRMATORIO.json` | `validation/p1/freeze_manifest.json` | Rastreabilidade e hashes | Selecionado; verificar referências/paths antes de copiar |
| 1 | `05_Experimentos_LLM/CAMPANHA_CONFIRMATORIA_20261019/RESULTADOS_CIENTIFICOS_CONFIRMATORIA.md` | `docs/p1_results.md` | Interpretação dos resultados | Selecionado; revisão científica antes da abertura pública |
| 2 | `05_Experimentos_LLM/CAMPANHA_CONFIRMATORIA_20261019/MATRIZ_CONFIRMATORIA.csv` | `experiments/p1/run_matrix.csv` | Ordem e pareamento | Candidato; inspecionar |
| 2 | `05_Experimentos_LLM/CAMPANHA_CONFIRMATORIA_20261019/COMPARACAO_PARES_CONFIRMATORIA.csv` | `experiments/p1/paired_comparison.csv` | Evidência agregada | Candidato; inspecionar |
| 2 | `01_Obsidian/Artigo_SMNE/03_Modelo_Reduzido.md` | `docs/model_reduction.md` | Contexto físico e escolhas | Candidato; contém texto em revisão, não copiar integralmente |
| 2 | `01_Obsidian/Artigo_SMNE/05_Metodologia_Experimental.md` | `docs/methodology.md` | Metodologia | Candidato; mistura planejamento antigo e atualização |
| 3 | `01_Obsidian/Research/BENCHMARK_P2_ARTIGO.md` | `docs/p2_development_notes.md` | Desenvolvimento P2 | Reter na origem nesta fase |
| 3 | `03_Modelo_Reduzido/P2_ARTIGO_20261022/` | `benchmark/p2/` | Dados e análises P2 | Reter na origem até revisão do conjunto congelado |

## Regras para o piloto Open Math OS

1. Compartilhar inicialmente apenas a especificação P1 e os critérios públicos de avaliação, **sem fornecer solução de referência, resultados confirmatórios, respostas de LLM, códigos completos do avaliador ou testes privados**.
2. Dar ao auditor resultados gerados em outra etapa, **após** congelar a avaliação cega, se quisermos medir identificação de erros sem vazamento de respostas.
3. Preservar cópias imutáveis, SHA-256, proveniência, unidades, versão do solver e ambientes.
4. Não publicar nem modificar artefatos originais da campanha confirmatória; toda seleção é feita por cópia revisada.
5. Antes de alterar visibilidade do repositório, revisar direitos de redistribuição do deck Egg Model/Equinor, materiais de terceiros, licenças, dados e eventuais detalhes de ambiente.
6. Para uma auditoria independente, diferenciar **correção matemática**, **validade do surrogate** e **validação física no OPM Flow**.

## Fatos confirmados pela inspeção

- Protocolo P1: 8 pares, 16 sessões, A linguagem natural versus B SDD+LLM+TDD; quatro controles `INJECT3, INJECT6, INJECT7, INJECT8`, quatro controles fixos; T=3600 dias.
- O relatório de resultados P1 indica validade semântica 8/8 em cada condição, igualdade na correção matemática observada, documentação de rastreabilidade mais frequente em B e tempo maior em B. Essas afirmações vêm do relatório original, ainda não de reprodução independente neste repositório.
- P2 é um desenvolvimento separado: sete controles livres e surrogates de FOPT/FWIT; não importar decisões antigas como se fossem o protocolo final.

## Pendências antes de copiar

- Conferir o snapshot/commit exato do repositório de origem e a integridade dos arquivos contra o congelamento.
- Inspecionar dados derivados, prompts, rubrica pública e dependências mínimas para reprodução.
- Definir pacote **cego** para teste no Open Math OS versus pacote **público de reprodução** do artigo.
- Conferir o estado científico final do P2 e a versão atual do manuscrito.
