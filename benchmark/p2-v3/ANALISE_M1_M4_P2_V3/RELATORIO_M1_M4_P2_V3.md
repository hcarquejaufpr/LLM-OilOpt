# P2 v3 — análise M1–M4

Análise pós-unblinding determinística. A avaliação semântica cega foi reutilizada sem alterações.

| Métrica | A | B |
|---|---:|---:|
| M1 semantic PASS | 7 | 8 |
| M2 candidatos 5/5 execução | 7 | 8 |
| M2 média políticas executadas | 4.375 | 5.0 |
| M2 candidatos 5/5 ref. numérica | 7 | 8 |
| M3 mediana MAE FOPT | 6.111804395914078e-09 | 6.1583705246448515e-09 |
| M3 mediana MAE FWIT | 4.1727907955646516e-07 | 4.1727907955646516e-07 |
| M3 mediana max |erro q| | 8.161165965248074e-09 | 8.161165965248074e-09 |
| M3 mediana max |erro r| | 7.738421459180245e-09 | 7.738421459180245e-09 |
| M4 ranking PASS | 7 | 8 |
| M4 trade-off PASS | 7 | 8 |

## Pareamento

| Par | A | B | M1 A/B | M2 A/B | Δ(B−A) MAE FOPT | Δ(B−A) MAE FWIT |
|---|---|---|---|---:|---:|---:|
| CF01 | P2V3-C13 | P2V3-C12 | True/True | 5/5 | 1.711305230855942e-09 | 5.946494638919827e-08 |
| CF02 | P2V3-C01 | P2V3-C05 | False/True | 0/5 | NA | NA |
| CF03 | P2V3-C14 | P2V3-C07 | True/True | 5/5 | -2.156768459826708e-06 | -6.789867766201496e-05 |
| CF04 | P2V3-C09 | P2V3-C16 | True/True | 5/5 | -6.984919309616105e-11 | 1.3969838619235354e-10 |
| CF05 | P2V3-C06 | P2V3-C02 | True/True | 5/5 | 2.280576154589653e-08 | 6.094109266996385e-07 |
| CF06 | P2V3-C10 | P2V3-C08 | True/True | 5/5 | 5.9371814131736755e-09 | 4.1527673602104185e-07 |
| CF07 | P2V3-C04 | P2V3-C03 | True/True | 5/5 | -1.062871888279915e-08 | 9.313225746156903e-11 |
| CF08 | P2V3-C11 | P2V3-C15 | True/True | 5/5 | -9.313225746154835e-11 | 0.0 |

## Inferência pareada exploratória

- M1 McNemar exato: p=1.0
- M2 McNemar exato: p=1.0
- M3 FOPT MAE sign test: p=1.0
- M3 FWIT MAE sign test: p=0.21875

**Guardrail:** n=8 pares; interpretar como evidência exploratória, não como prova de superioridade.
