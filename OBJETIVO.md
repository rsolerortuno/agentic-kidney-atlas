# OBJETIVO: GraphRAG multiómico ccRCC

Construir y dejar funcionando, de principio a fin, un sistema GraphRAG multiómico para carcinoma renal de células claras (ccRCC). A partir de datos bulk de CPTAC + literatura (PubTator3) + APIs externas vía MCP, un orquestador y 3 subagentes de Claude Code deben responder "¿Cuáles son las dianas terapéuticas y biomarcadores candidatos en ccRCC con mayor evidencia convergente?" y generar un informe HTML trazable. Single-cell y espacial quedan FUERA de este objetivo.

## Restricciones
- Windows 11, sin Python/Docker/Java. Autorizado: `winget install astral-sh.uv` (solo usuario). No instalar nada más a nivel de sistema.
- Proyecto Python 3.12 con uv, en git (un commit por milestone).
- LLM: solo Claude Code con la suscripción Pro. Nada de pago aparte: no usar la API de Claude, el SDK `anthropic`, Claude Agent SDK ni ANTHROPIC_API_KEY. Todo el razonamiento LLM lo hacen Claude Code y sus subagentes:
  - Orquestador: sesión principal, Claude Opus 5.5, esfuerzo medium.
  - Todos los subagentes (agente-omico, agente-literatura, agente-externo, extractor-relaciones): `model: sonnet` (Sonnet 5.5).
- Ahorro de cuota: extracción de relaciones limitada a 100 abstracts; el resto de la literatura, solo con anotaciones de PubTator3. Nunca pasar matrices o el grafo entero en un prompt: los agentes consultan por herramientas MCP.
- Grafo embebido: SQLite + NetworkX detrás de una interfaz `GraphStore`, para migrar a Neo4j sin tocar agentes ni informe.
- Todo cálculo numérico (DE, puntuación de convergencia, verificación) es código Python determinista, nunca LLM.
- Ningún hallazgo sin evidence_id. Nada de valores inventados. Si algo no se puede verificar, se dice explícitamente.
- Estado del progreso en PROGRESS.md (milestone actual, comprobaciones superadas, bloqueos).

## Milestones (cada uno se cierra solo cuando su comprobación pasa)

**M0 Entorno:** uv + pyproject (pandas, pyarrow, duckdb, scipy, statsmodels, networkx, requests, pydantic, mcp/fastmcp, jinja2, plotly, rank_bm25, pytest).
- ✔ `uv run python -c "import pandas, mcp, networkx"` funciona; `uv run pytest` se ejecuta.

**M1 Datos:** descargar de https://linkedomics.org/data_download/CPTAC-CCRCC/ los 8 archivos (HS_CPTAC_CCRCC_CLI.tsi, HS_CPTAC_CCRCC_CNV_gene_Tumor.cct, HS_CPTAC_CCRCC_proteome_Tumor.cct, HS_CPTAC_CCRCC_proteome_Normal.cct, HS_CPTAC_CCRCC_phosphoproteome_gene_Tumor.cct, HS_CPTAC_CCRCC_phosphoproteome_gene_Normal.cct, HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct, HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct) a data/raw/ + conjunto completo de HGNC. Manifest con URL, fecha, sha256 y dimensiones.
- ✔ ≈110 tumores; 84 normales en proteoma y fosfoproteoma, 75 en RNA.

**M2 Análisis ómico:** mapear a HGNC/Ensembl; DE tumor frente a normal (emparejado por paciente cuando sea posible; t pareado o Wilcoxon + Benjamini-Hochberg) en RNA, proteoma y fosfoproteoma; frecuencia de ganancias y pérdidas en CNV; correlación RNA–proteína por gen. Salida en data/processed/*.parquet.
- ✔ CA9, NDUFA4L2, ANGPTL4, EGLN3 y VEGFA significativamente aumentados (padj<0.05) en RNA y en proteína cuando están medidos; pérdida de 3p frecuente. tests/test_omics.py en verde.

**M3 Literatura:** PubTator3 API → PMIDs de ccRCC + top 50 genes de M2, con anotaciones normalizadas. El subagente extractor-relaciones procesa como máximo 100 abstracts con esquema cerrado {sujeto, relación ∈ [upregulated_in, downregulated_in, biomarker_of, target_of, drives, associated_with, inhibits], objeto, frase_literal, pmid, confianza} y escribe JSON a data/literature/. Índice BM25 de fragmentos.
- ✔ Script que comprueba que cada frase_literal aparece en su abstract; se descartan las que no; 0 relaciones sin PMID.

**M4 Grafo de conocimiento:** nodos Gene, Protein, Pathway, Drug, Disease (MONDO:0005005), Dataset, Publication, Chunk, Evidence. Toda relación pasa por un nodo Evidence (tipo, estadístico, fuente, fecha). Cargar M2, M3 y la caché de M5.
- ✔ Informe de recuentos; el 100 % de las relaciones tiene ≥1 Evidence con fuente resoluble.

**M5 Servidores MCP** (FastMCP, stdio), registrados en .mcp.json:
- omics-mcp: get_de_results(gene, layer), top_genes(layer, n), rna_protein_correlation(gene).
- kg-mcp: neighbors(entity), evidence_for(entity), search_literature(query), check_evidence_ids(ids).
- external-mcp: Open Targets GraphQL (asociación diana–enfermedad, tractabilidad, fármacos), STRING (interacciones), UniProt (función). Respuestas cacheadas en disco.
- ✔ Un test pytest por herramienta, con llamadas reales y caché.

**M6 Agentes (Claude Code):**
- .claude/agents/agente-omico.md (model: sonnet): solo omics-mcp + kg-mcp.
- .claude/agents/agente-literatura.md (model: sonnet): solo kg-mcp.
- .claude/agents/agente-externo.md (model: sonnet): solo external-mcp.
- .claude/agents/extractor-relaciones.md (model: sonnet): usado en M3.
- .claude/commands/informe.md (orquestador, sesión principal con Opus 5.5 medium): lanza los 3 agentes en paralelo con la pregunta; cada uno escribe JSON con evidence_ids en runs/<fecha>/; luego ejecuta scripts/integrate.py (puntuación de convergencia determinista, pesos en config.yaml) y scripts/verify.py (descarta y registra los evidence_ids inexistentes); por último redacta los párrafos de interpretación solo a partir del JSON verificado y escribe runs/<fecha>/findings.json.
- ✔ `claude -p "/informe" --model claude-opus-5-5` genera findings.json; verify.py reporta 0 afirmaciones sin soporte.

**M7 Informe HTML:** plantilla Jinja2 determinista (scripts/render.py) a partir de findings.json. Secciones: resumen, hallazgos por capa (volcano en Plotly), tabla de convergencia (top 25), red interactiva (Cytoscape.js), dianas con fármacos (Open Targets), citas enlazadas a PubMed, procedencia (datasets, versiones, fecha, modelos usados) y limitaciones.
- ✔ Un único reports/ccrcc_report.html que se abre sin servidor; un validador comprueba que cada afirmación enlaza a un evidence_id existente.

**M8 Evaluación:** (a) recall de VHL, PBRM1, SETD2, BAP1, KDM5C, EPAS1, CA9, VEGFA, MTOR, HIF1A + top 20 de Open Targets para MONDO:0005005, medido en el top 100 de convergencia; (b) 10 preguntas con respuesta conocida en evals/qa.yaml, respondidas con `claude -p`.
- ✔ Recall ≥ 60 % y ≥ 8/10 correctas, o un análisis de por qué no y qué se cambió. Resultados en evals/results.md.

**M9 Reproducibilidad:** `uv run ccrcc run-all` ejecuta M1→M5 y el render; README con arquitectura, uso, consumo aproximado de cuota, cómo pasar a Neo4j y cómo añadir sc/espacial (Li 2022 Mendeley g67bkbnhhg, GSE175540, GSE156632).
- ✔ Una ejecución limpia regenera el informe sin errores.

## Definición de terminado
M0–M9 superados, `uv run pytest` en verde, reports/ccrcc_report.html generado con datos reales y evaluado en evals/results.md. Si se agota la cuota de Pro a mitad, guardar el estado en PROGRESS.md, anotar en qué milestone se quedó y continuar desde ahí en la siguiente ventana de uso.
