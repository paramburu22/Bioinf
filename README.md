# Bioinf - Trabajo Práctico

Repositorio con los scripts, datos y documentación del TP de Bioinformática enfocado en el gen **HBB**. Todo el material fue normalizado para facilitar la ejecución de los ejercicios y la revisión del código.

---

## 📁 Estructura

```
Bioinf/
├── config/                     # Archivos de configuración (JSON, etc.)
├── data/
│   ├── raw/                    # Datos fuente (GenBank, etc.)
│   ├── interim/                # Productos intermedios (ORFs, CDS, ...)
│   ├── results/
│   │   ├── ex02/               # Resultados del BLAST remoto
│   │   ├── ex03/               # Insumos/resultados del MSA
│   │   ├── ex04/               # Reportes EMBOSS / PROSITE
│   │   └── ex05/               # Primers generados
│   └── external/prosite/       # Bases de datos descargadas
├── docs/                       # Material para la exposición del TP
├── scripts/                    # Wrappers Bash para cada ejercicio
├── src/
│   ├── exercises/              # Scripts Python principales (ex01..ex05)
│   └── utils/                  # Utilidades y ayudas
└── README.md
```

---

## ⚙️ Configuración Inicial

1. **Dependencias Python**
   ```bash
   pip3 install -r requirements.txt
   ```

2. **EMBOSS (para Ej. 4)**
   - Usar Conda / Bioconda o el helper `scripts/install_emboss.sh` si trabajas en macOS.

---

## 🧪 Ejercicios

### Ejercicio 1 – ORFs desde GenBank
- Script: `src/exercises/ex01_generate_orfs.py`
- Input: `data/raw/hbb_nm000518.gb`
- Output: `data/interim/hbb_orfs.fasta`
- Ejecución:
  ```bash
  python3 src/exercises/ex01_generate_orfs.py
  ```

### Ejercicio 2 – BLAST remoto + parser
- Script BLAST: `src/exercises/ex02_blast_remote.py`
- Parser opcional: `src/exercises/ex02_blast_parser.py`
- Input: `data/interim/hbb_orfs.fasta`
- Output principal: `data/results/ex02/blast_hbb_remote.xml`
- Ejemplo:
  ```bash
  python3 src/exercises/ex02_blast_remote.py --hitlist-size 15
  python3 src/exercises/ex02_blast_parser.py --results-dir data/results/ex02
  ```

### Ejercicio 3 – MSA con Biopython
- Script Bash: `scripts/run_ex3.sh`
- Script Python: `src/exercises/ex03_msa.py`
- Inputs:
  - `data/interim/hbb_correct_orf.fasta` (generado con `src/utils/extract_correct_orf.py`)
  - `data/results/ex02/blast_hbb_remote.xml`
- Salidas: `data/results/ex03/{msa_input.fasta, msa_alignment.fasta, msa_summary.txt}`
- Ejecución:
  ```bash
  # ENTREZ_EMAIL es obligatorio (argumento o variable de entorno)
  ./scripts/run_ex3.sh tu_email@institucion.edu --top 12
  ```

### Ejercicio 4 – PROSITE / EMBOSS
- Script Bash: `scripts/run_ex4.sh`
- Script Python: `src/exercises/ex04_emboss_prosite.py`
- Inputs: `data/interim/hbb_orfs.fasta` + base PROSITE en `data/external/prosite`
- Output: `data/results/ex04/hbb_domain_analysis.txt`
- Ejecuta:
  ```bash
  ./scripts/run_ex4.sh
  ```

### Ejercicio 5 – Diseño de primers
- Script Bash: `scripts/run_ex5.sh`
- Script Python: `src/exercises/ex05_primer_design.py`
- Inputs:
  - `data/raw/hbb_nm000518.gb`
  - `config/primer_design.json` (se crea automáticamente si no existe)
- Output: `data/results/ex05/hbb_primers.fasta`
- Manual:
  ```bash
  python3 src/exercises/ex05_primer_design.py
  ```

---

## 📂 Datos Clave

| Archivo/Directorio                             | Descripción                                |
|------------------------------------------------|--------------------------------------------|
| `data/raw/hbb_nm000518.gb`                     | GenBank del transcript HBB (NM_000518)     |
| `data/interim/hbb_orfs.fasta`                  | ORFs de los seis marcos (Ej.1)             |
| `data/interim/hbb_correct_orf.fasta`           | CDS traducida (marco correcto)             |
| `data/results/ex02/blast_hbb_remote.xml`       | BLASTP vs SwissProt (Ej.2)                 |
| `data/results/ex03/msa_summary.txt`            | Informe del alineamiento múltiple          |
| `data/results/ex04/hbb_domain_analysis.txt`    | Reporte PROSITE / EMBOSS                   |
| `data/results/ex05/hbb_primers.fasta`          | FASTA con los primers seleccionados        |

---

## 📚 Documentación Complementaria

- `docs/exercise03_overview.md`: resumen del pipeline MSA + hallazgos.
- `docs/exercise04_overview.md`: procedimiento y resultados EMBOSS.
- `docs/exercise05_overview.md`: criterios de diseño de primers.
- `docs/exercise06_research.md`: notas de la investigación bibliográfica.

---

## 🗒️ Notas finales

- Todos los scripts usan rutas relativas al **root del repositorio**, por lo que conviene ejecutar los comandos desde `Bioinf/` (los wrappers Bash ya lo hacen).
- Para reutilizar el pipeline con otro gen, basta con reemplazar los datos en `data/raw` y ajustar los parámetros de cada script/JSON.
- Cualquier salida adicional debería guardarse dentro de `data/results/<ejercicio>` para mantener la estructura ordenada.
