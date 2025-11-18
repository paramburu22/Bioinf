## Ejercicio 3 – Alineamiento múltiple (HBB)

### Flujo implementado
1. **Entrada BLAST:** se reutiliza el resultado remoto generado en el Ej.2 (`data/results/ex02/blast_hbb_remote.xml`) para identificar los 10 mejores hits contra SwissProt.
2. **Descarga automatizada:** el script `src/exercises/ex03_msa.py` (envuelto por `scripts/run_ex3.sh`) usa `Bio.Entrez` para bajar la secuencia proteica de cada accession seleccionado. Es obligatorio brindar un email válido (argumento del `.sh` o variable `ENTREZ_EMAIL`).
3. **Construcción del FASTA:** se concatena la secuencia correcta del CDS (`data/interim/hbb_correct_orf.fasta`) con los hits descargados → `data/results/ex03/msa_input.fasta`.
4. **MSA con Biopython:** se generan alineamientos globales pareados (`Bio.pairwise2`) contra la secuencia consulta y se combinan en un star-MSA sin depender de binarios externos. El resultado se guarda en `data/results/ex03/msa_alignment.fasta`.
5. **Interpretación:** se calcula automáticamente una tabla de identidades, el porcentaje de columnas conservadas y un consenso con umbral 70%. El informe se guarda en `data/results/ex03/msa_summary.txt`.

### Cómo ejecutarlo

```bash
# Instalar dependencias (una vez)
pip3 install -r requirements.txt

# Ejecutar el Ejercicio 3
./scripts/run_ex3.sh tu_email@institucion.edu
```

Parámetros adicionales (por ejemplo usar otro XML):

```bash
./scripts/run_ex3.sh tu_email@institucion.edu --blast-xml data/results/ex02/otro.xml --top 8
```

### Resultados destacados
- **Identidad promedio** frente a la beta-globina humana: `~97%` entre los 10 mejores hits.
- **Conservación**: 135 de 147 posiciones están 100% conservadas; 140 posiciones ≥80%.
- **Consenso**: retiene las histidinas proximal/distal implicadas en la unión del grupo hemo y evidencia la alta conservación del pliegue globina.

Los archivos producidos (`data/results/ex03/msa_input.fasta`, `data/results/ex03/msa_alignment.fasta`, `data/results/ex03/msa_summary.txt`) sirven como insumo para la discusión solicitada en la consigna.

