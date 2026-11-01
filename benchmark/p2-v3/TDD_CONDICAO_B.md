# TDD público — Condição B (P2 v3)

Use `B/public_tdd.py` como suíte pública mínima. Ela verifica apenas contrato e invariantes visíveis e não contém o ótimo nem resultados OPM.

Ciclo esperado:
1. execute os testes antes da implementação e registre a falha quando aplicável;
2. implemente os requisitos;
3. execute novamente até obter green;
4. registre em `evidence.json` as evidências R01–R12.

Além do contrato básico, o teste público verifica que FOPT reportado é coerente com a definição pública v3, inclusive `real_i=(r_i-80)/20`.

Os testes públicos não substituem o evaluator privado.
