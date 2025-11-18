# Ejercicio 4 - Análisis de Dominios PROSITE con EMBOSS

## Descripción General

El Ejercicio 4 consiste en realizar un análisis de dominios de proteínas utilizando EMBOSS y la base de datos PROSITE. El objetivo es identificar motivos funcionales en las secuencias de aminoácidos obtenidas del Ejercicio 1.

## ¿Qué es PROSITE?

PROSITE es una base de datos de patrones y perfiles de proteínas que permite identificar:
- **Dominios funcionales**: regiones de proteínas con funciones específicas
- **Sitios de modificación post-traduccional**: lugares donde las proteínas pueden ser modificadas
- **Motivos estructurales**: patrones de secuencia conservados

## ¿Qué es EMBOSS?

EMBOSS (European Molecular Biology Open Software Suite) es un conjunto de herramientas de análisis bioinformático. Para este ejercicio utilizamos:
- **`patmatmotifs`**: Busca patrones PROSITE en secuencias de proteínas
- **`prosextract`**: Prepara la base de datos PROSITE para su uso

## Flujo del Ejercicio 4

### 1. Entrada: Archivo FASTA de Proteínas
- **Archivo**: `data/interim/hbb_orfs.fasta` (generado en el Ejercicio 1)
- **Contenido**: 6 secuencias de proteínas (6 marcos de lectura)
- **Origen**: Traducción de las 6 lecturas posibles del transcript HBB

### 2. Proceso de Análisis

#### Paso 1: Descarga de Base de Datos PROSITE
- Se descarga `data/external/prosite/prosite.dat` (patrones PROSITE)
- Se descarga `data/external/prosite/prosite.doc` (documentación, necesaria para `prosextract`)

#### Paso 2: Preparación de la Base de Datos
- Se ejecuta `prosextract` para preparar la base de datos PROSITE
- Esto permite que `patmatmotifs` pueda buscar patrones eficientemente

#### Paso 3: Análisis de Dominios
- Para cada secuencia de proteína en el archivo FASTA:
  - Se ejecuta `patmatmotifs` con las siguientes opciones:
    - `-full`: Obtener documentación completa de los motivos
    - `-noprune`: Incluir patrones simples (no ignorar modificaciones post-traduccionales)
  - Se buscan todos los motivos PROSITE que coinciden con la secuencia

### 3. Salida: Archivo de Resultados
- **Archivo**: `data/results/ex04/hbb_domain_analysis.txt`
- **Formato**: Texto con análisis detallado de cada secuencia

## Explicación de los Resultados

### Estructura del Archivo de Resultados

El archivo `data/results/ex04/hbb_domain_analysis.txt` contiene:

1. **Encabezado**: Información general del análisis
   - Archivo de entrada
   - Número de secuencias analizadas

2. **Análisis por Secuencia**: Para cada una de las 6 secuencias:
   - Nombre de la secuencia (ej: `NM_000518.5_frame1`)
   - Longitud en aminoácidos
   - Lista de motivos PROSITE encontrados

### Tipos de Motivos PROSITE Encontrados

#### 1. Sitios de Fosforilación

**PKC_PHOSPHO_SITE**
- Sitio de fosforilación por Protein Kinase C
- Patrón: `[ST]-x-[RK]`
- Función: Regulación de actividad proteica mediante fosforilación

**CK2_PHOSPHO_SITE**
- Sitio de fosforilación por Casein Kinase 2
- Patrón: `[ST]-x(2)-[DE]`
- Función: Modificación post-traduccional común

**CAMP_PHOSPHO_SITE**
- Sitio de fosforilación dependiente de cAMP
- Patrón: `[RK](2)-x-[ST]`
- Función: Regulación por vía de señalización cAMP

#### 2. Modificaciones Post-Traduccionales

**MYRISTYL**
- Sitio de miristoilación
- Patrón: `G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}`
- Función: Unión de ácido mirístico (anclaje a membrana)

**ASN_GLYCOSYLATION**
- Sitio de glicosilación en asparagina
- Patrón: `N-{P}-[ST]-{P}`
- Función: Adición de carbohidratos a la proteína

**AMIDATION**
- Sitio de amidación
- Patrón: `[GKRHDENQSTAVILM]-x-[KRHDENQSTAVILM]-G-[KRHDENQSTAVILM]-[KR]-[KRHDENQSTAVILM]`
- Función: Procesamiento de péptidos

### Interpretación de los Resultados

#### Resumen por Frame:

| Frame | Motivos Encontrados | Observaciones |
|-------|---------------------|---------------|
| Frame 1 | 6 motivos | Incluye fosforilación (PKC, cAMP) y miristoilación |
| Frame 2 | 10 motivos | Mayor diversidad de motivos |
| Frame 3 | 8 motivos | **Incluye AMIDATION** (relevante biológicamente) |
| Frame 4 | 13 motivos | Mayor cantidad; incluye glicosilación |
| Frame 5 | 7 motivos | Varios tipos de fosforilación |
| Frame 6 | 6 motivos | Sitios de fosforilación |

#### ¿Por qué Frame 3 es el más relevante?

- **Frame 3** corresponde probablemente a la **proteína funcional real** (hemoglobina beta)
- Contiene el motivo **AMIDATION**, que es biológicamente relevante para hemoglobina
- Los otros frames son artefactos de las traducciones alternativas, pero el análisis PROSITE los detecta igualmente

### Información Detallada de Cada Motivo

Para cada motivo encontrado, el archivo muestra:

1. **Nombre del Motivo**: Tipo de dominio/motivo PROSITE
2. **Secuencia del Contexto**: Secuencia de aminoácidos alrededor del sitio
3. **Longitud**: Tamaño del motivo en aminoácidos
4. **Posiciones**: Start y End en la secuencia

**Ejemplo:**
```
Motivo: AMIDATION
PKVKAHGKKVLGAF
     |  |
  80  83
```
- El motivo AMIDATION está en las posiciones 80-83
- La secuencia muestra el contexto alrededor del sitio
- Los `|` indican el sitio exacto del motivo

## Uso del Script

### Ejecución

```bash
# Opción 1: Usar el script de instalación y ejecución
./scripts/run_ex4.sh

# Opción 2: Ejecutar directamente
python ex4_emboss_analysis.py
```

### Requisitos

- **Python 3.x** con Biopython
- **EMBOSS** instalado (se puede instalar con conda: `conda install -c bioconda emboss`)
- **Archivo de entrada**: `data/interim/hbb_orfs.fasta` (generado por el Ejercicio 1)

### Archivos Generados

1. **`data/results/ex04/hbb_domain_analysis.txt`**: Resultados del análisis de dominios
2. **`prosite.dat`**: Base de datos PROSITE (descargada automáticamente)
3. **`prosite.doc`**: Documentación PROSITE (descargada automáticamente)
4. **`prosite/`**: Directorio con base de datos preparada por `prosextract`

## Importancia Biológica

### ¿Por qué es importante analizar dominios PROSITE?

1. **Predicción de Función**: Los dominios PROSITE pueden indicar la función de una proteína
2. **Modificaciones Post-Traduccionales**: Identifica sitios donde la proteína puede ser modificada
3. **Regulación**: Los sitios de fosforilación indican posibles puntos de regulación
4. **Análisis Comparativo**: Permite comparar proteínas y encontrar similitudes funcionales

### En el Contexto de HBB (Hemoglobina Beta)

- **AMIDATION**: Procesamiento de péptidos, relevante para la maduración de hemoglobina
- **Sitios de Fosforilación**: Pueden regular la función de la hemoglobina
- **Glicosilación**: Puede afectar la estabilidad y función de la proteína

## Limitaciones y Consideraciones

1. **Falsos Positivos**: No todos los motivos detectados son biológicamente relevantes
2. **Frames Alternativos**: Los frames 1, 2, 4, 5 y 6 son artefactos de traducción
3. **Patrones Simples**: Algunos motivos son muy simples y pueden aparecer por casualidad
4. **Validación**: Los resultados deben validarse con análisis experimentales

## Conclusión

El Ejercicio 4 demuestra cómo:
- Utilizar herramientas bioinformáticas (EMBOSS) para análisis de proteínas
- Identificar dominios funcionales usando bases de datos especializadas (PROSITE)
- Interpretar resultados de análisis de secuencias de proteínas
- Relacionar patrones de secuencia con función biológica

El archivo `data/results/ex04/hbb_domain_analysis.txt` contiene un análisis completo de dominios PROSITE en las secuencias de proteínas del Ejercicio 1, identificando 50 motivos en total distribuidos en las 6 secuencias analizadas.

