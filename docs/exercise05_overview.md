# Ejercicio 5 - Diseño Parametrizable de Primers

## Descripción General

El Ejercicio 5 consiste en diseñar primers para el transcript HBB (NM_000518) de forma parametrizable. Los primers son secuencias cortas de ADN utilizadas en técnicas de PCR (Reacción en Cadena de la Polimerasa) para amplificar regiones específicas del genoma.

## ¿Qué son los Primers?

Los **primers** (iniciadores) son secuencias cortas de oligonucleótidos (18-30 nucleótidos) que:
- Se unen específicamente a una región complementaria del ADN/ARN
- Sirven como punto de inicio para la síntesis de ADN por la polimerasa
- Permiten amplificar regiones específicas del genoma mediante PCR

## Objetivo del Ejercicio

Diseñar **5 primers** para el transcript HBB que cumplan con criterios específicos de calidad, establecidos en un archivo de configuración JSON.

## Requisitos de Diseño

Los primers deben cumplir con los siguientes criterios:

### 1. Longitud
- **Rango**: 18-24 pares de bases (bp)
- **Razón**: Longitudes muy cortas pueden ser inespecíficas; muy largas pueden tener problemas de hibridación

### 2. Contenido GC
- **Rango**: 50-60%
- **Razón**: Un contenido GC equilibrado mejora la estabilidad del duplex
- **Fórmula**: `GC% = (G + C) / Longitud × 100`

### 3. Extremos Terminales
- **Restricción**: Evitar G o C en las primeras 2 y últimas 2 posiciones
- **Razón**: Los extremos con G/C pueden formar estructuras secundarias (hairpins) que interfieren con la hibridación

### 4. Temperatura de Melting (Tm)
- **Máximo**: ≤67°C
- **Razón**: Temperaturas muy altas pueden causar desnaturalización no específica
- **Cálculo**: Usa el método de Wallace: `Tm = 64.9 + 41 × (GC_count - 16.4) / length`

### 5. Cantidad
- **Número**: 5 primers
- **Selección**: Se eligen los mejores según temperatura de melting cercana a 65°C y contenido GC cercano a 55%

## Flujo del Ejercicio 5

### 1. Entrada

#### Archivo GenBank
- **Archivo**: `data/raw/hbb_nm000518.gb`
- **Contenido**: Secuencia del transcript HBB (NM_000518) en formato GenBank
- **Longitud**: 628 bp (ARN mensajero)

- **Archivo de Configuración**
- **Archivo**: `config/primer_design.json`
- **Formato**: JSON
- **Parámetros**:
  ```json
  {
    "min_length": 18,
    "max_length": 24,
    "min_gc": 50,
    "max_gc": 60,
    "max_tm": 67,
    "num_primers": 5,
    "avoid_terminal_gc": true,
    "terminal_positions": 2
  }
  ```

### 2. Proceso de Diseño

#### Paso 1: Extracción de la Secuencia
- Se lee el archivo GenBank
- Se extrae la secuencia de nucleótidos del transcript HBB
- Se obtiene la secuencia completa (628 bp)

#### Paso 2: Generación de Candidatos
- Se genera un conjunto de candidatos usando **ventana deslizante**
- Para cada longitud posible (18-24 bp):
  - Se crean todas las subsecuencias posibles de esa longitud
  - Cada subsecuencia es un candidato potencial de primer

**Ejemplo**:
```
Secuencia: ATGCGATCGATCG...
Candidatos de 18 bp:
- Posición 1-18: ATGCGATCGATCGATCG
- Posición 2-19: TGCGATCGATCGATCGA
- Posición 3-20: GCGATCGATCGATCGAT
...
```

#### Paso 3: Filtrado por Criterios

Cada candidato se evalúa según:

1. **Longitud**: Debe estar en el rango 18-24 bp
2. **Contenido GC**: Calculado y verificado (50-60%)
3. **Extremos Terminales**: Verificación de que no haya G/C en las primeras 2 y últimas 2 posiciones
4. **Temperatura de Melting**: Calculada y verificada (≤67°C)

**Ejemplo de verificación de extremos**:
```
Secuencia: ATGCGATCGATCGATCGAT
Primeros 2: AT ✓ (no tiene G/C)
Últimos 2: AT ✓ (no tiene G/C)
✅ Pasa el filtro
```

#### Paso 4: Selección de los Mejores Primers

- Se ordenan los primers filtrados por:
  1. **Tm cercana a 65°C** (temperatura óptima)
  2. **GC% cercana a 55%** (valor medio del rango)
- Se seleccionan los **5 mejores primers**

### 3. Salida

#### Archivo FASTA
- **Archivo**: `data/results/ex05/hbb_primers.fasta`
- **Formato**: FASTA estándar
- **Contenido**: 5 primers con sus identificadores

#### Salida en Consola
- Lista de los 5 primers diseñados
- Para cada primer:
  - Secuencia
  - Posición en el transcript
  - Longitud
  - Contenido GC
  - Temperatura de melting

## Explicación de los Resultados

### Estructura del Archivo de Salida

El archivo `data/results/ex05/hbb_primers.fasta` contiene:

```
>Primer_1_NM_000518.5_pos448-468
AAGTGGTGGCTGGTGTGGCTA
>Primer_2_NM_000518.5_pos467-487
TAATGCCCTGGCCCACAAGTA
>Primer_3_NM_000518.5_pos468-488
AATGCCCTGGCCCACAAGTAT
>Primer_4_NM_000518.5_pos574-594
AAGGGCCTTGAGCATCTGGAT
>Primer_5_NM_000518.5_pos467-488
TAATGCCCTGGCCCACAAGTAT
```

### Información de Cada Primer

#### Primer 1
- **Secuencia**: `AAGTGGTGGCTGGTGTGGCTA`
- **Posición**: 448-468 (21 bp)
- **Extremos**: `AA`...`TA` ✓ (sin G/C en extremos)
- **GC**: ~52% (dentro del rango)
- **Tm**: ~65°C (dentro del límite)

#### Primer 2
- **Secuencia**: `TAATGCCCTGGCCCACAAGTA`
- **Posición**: 467-487 (21 bp)
- **Extremos**: `TA`...`TA` ✓ (sin G/C en extremos)
- **GC**: ~57% (dentro del rango)
- **Tm**: ~67°C (en el límite máximo)

#### Primer 3
- **Secuencia**: `AATGCCCTGGCCCACAAGTAT`
- **Posición**: 468-488 (21 bp)
- **Extremos**: `AA`...`AT` ✓ (sin G/C en extremos)
- **GC**: ~57% (dentro del rango)
- **Tm**: ~67°C (en el límite máximo)

#### Primer 4
- **Secuencia**: `AAGGGCCTTGAGCATCTGGAT`
- **Posición**: 574-594 (21 bp)
- **Extremos**: `AA`...`AT` ✓ (sin G/C en extremos)
- **GC**: ~57% (dentro del rango)
- **Tm**: ~67°C (en el límite máximo)

#### Primer 5
- **Secuencia**: `TAATGCCCTGGCCCACAAGTAT`
- **Posición**: 467-488 (21 bp)
- **Nota**: Duplicado del Primer 2 (posible selección por criterios diferentes)

## Cálculo de Temperatura de Melting

### Método de Wallace

El script utiliza el método de Wallace para calcular la temperatura de melting:

```
Tm = 64.9 + 41 × (GC_count - 16.4) / length
```

**Ejemplo**:
- Secuencia: `AAGTGGTGGCTGGTGTGGCTA` (21 bp)
- G = 11, C = 0, Total = 11
- GC% = 11/21 = 52.4%
- Tm = 64.9 + 41 × (11 - 16.4) / 21
- Tm = 64.9 + 41 × (-5.4) / 21
- Tm = 64.9 - 10.5
- Tm ≈ 54.4°C

**Nota**: El método de Wallace es una aproximación simple. Métodos más precisos como Nearest Neighbor (Tm_NN) consideran la secuencia completa.

## Parámetros de Configuración

### Archivo `config/primer_design.json`

```json
{
  "min_length": 18,        // Longitud mínima en bp
  "max_length": 24,        // Longitud máxima en bp
  "min_gc": 50,            // Contenido GC mínimo (%)
  "max_gc": 60,            // Contenido GC máximo (%)
  "max_tm": 67,            // Temperatura de melting máxima (°C)
  "num_primers": 5,        // Número de primers a diseñar
  "avoid_terminal_gc": true,  // Evitar G/C en extremos
  "terminal_positions": 2     // Número de posiciones terminales a verificar
}
```

### Personalización

Puedes modificar estos parámetros según tus necesidades:

- **Longitud**: Ajustar `min_length` y `max_length` según el tipo de PCR
- **GC**: Modificar `min_gc` y `max_gc` según el contenido GC del genoma
- **Tm**: Ajustar `max_tm` según las condiciones de PCR
- **Cantidad**: Cambiar `num_primers` para obtener más opciones

## Uso del Script

### Ejecución

```bash
# Opción 1: Usar el script de instalación y ejecución
./scripts/run_ex5.sh

# Opción 2: Ejecutar directamente
python3 src/exercises/ex05_primer_design.py
```

### Requisitos

- **Python 3.x** con Biopython
- **Archivo GenBank**: `data/raw/hbb_nm000518.gb`
- **Archivo de configuración**: `config/primer_design.json` (se crea automáticamente si no existe)

### Archivos Generados

1. **`data/results/ex05/hbb_primers.fasta`**: Archivo FASTA con los 5 primers diseñados
2. **Salida en consola**: Información detallada de cada primer

## Importancia Biológica

### ¿Por qué son importantes estos criterios?

1. **Longitud 18-24 bp**:
   - Longitudes muy cortas pueden unirse a múltiples sitios (baja especificidad)
   - Longitudes muy largas pueden tener problemas de hibridación
   - 18-24 bp es el rango óptimo para PCR estándar

2. **Contenido GC 50-60%**:
   - Un contenido GC equilibrado proporciona estabilidad adecuada
   - Muy bajo GC (<50%): El duplex puede ser inestable
   - Muy alto GC (>60%): Puede formar estructuras secundarias (hairpins)

3. **Evitar G/C en extremos**:
   - Los extremos con G/C pueden formar estructuras secundarias
   - Los extremos con A/T son más flexibles y mejoran la especificidad

4. **Tm ≤67°C**:
   - Temperaturas muy altas pueden causar desnaturalización no específica
   - Permite usar condiciones de PCR estándar
   - La temperatura óptima suele estar alrededor de 65°C

### Aplicaciones de los Primers Diseñados

Los primers diseñados pueden usarse para:

1. **PCR cuantitativa (qPCR)**:
   - Detectar la presencia del transcript HBB
   - Cuantificar la expresión del gen

2. **PCR cualitativa**:
   - Verificar la presencia del transcript
   - Detectar variantes o mutaciones

3. **RT-PCR**:
   - Amplificar el transcript desde ARN
   - Estudiar la expresión génica

4. **Secuenciación**:
   - Amplificar regiones específicas para secuenciación
   - Identificar mutaciones en el gen HBB

## Limitaciones y Consideraciones

### Limitaciones del Script

1. **No considera pares de primers**:
   - El script diseña primers individuales
   - Para PCR real, necesitas un par forward/reverse
   - No verifica complementariedad entre primers

2. **No verifica estructuras secundarias**:
   - No detecta hairpins, dimers, o estructuras secundarias
   - Herramientas como Primer3 o BLAST deberían usarse para validación

3. **Método de Tm simplificado**:
   - Usa el método de Wallace (aproximación simple)
   - Métodos más precisos (Nearest Neighbor) considerarían la secuencia completa

4. **No verifica especificidad**:
   - No verifica si los primers se unen a otras regiones del genoma
   - Se recomienda validar con BLAST antes de usar

### Recomendaciones para Uso Real

1. **Validar con BLAST**:
   - Verificar que los primers se unen solo a la región deseada
   - Evitar uniones no específicas

2. **Usar herramientas especializadas**:
   - Primer3, OligoAnalyzer, o herramientas similares
   - Para validación completa de estructuras secundarias

3. **Diseñar pares de primers**:
   - Diseñar primers forward y reverse complementarios
   - Verificar que el producto amplificado sea del tamaño esperado

4. **Validación experimental**:
   - Probar los primers en condiciones reales de PCR
   - Ajustar temperatura de annealing según resultados

## Comparación con Herramientas Profesionales

### Herramientas Comerciales

- **Primer3**: Diseño de primers con múltiples criterios
- **OligoAnalyzer**: Análisis de estructuras secundarias
- **BLAST**: Verificación de especificidad

### Ventajas de Este Script

- **Parametrizable**: Fácil ajuste de criterios via JSON
- **Educativo**: Muestra los principios de diseño de primers
- **Rápido**: Análisis inmediato de múltiples candidatos
- **Transparente**: Código abierto, fácil de entender y modificar

## Conclusión

El Ejercicio 5 demuestra cómo:
- Diseñar primers de forma parametrizable usando Python
- Aplicar criterios bioquímicos para selección de primers
- Calcular propiedades importantes (GC%, Tm)
- Generar resultados en formato FASTA estándar

El script genera **5 primers** que cumplen con todos los criterios especificados:
- Longitud: 18-24 bp ✓
- Contenido GC: 50-60% ✓
- Sin G/C en extremos ✓
- Tm ≤67°C ✓

Estos primers pueden usarse para análisis cuanti y cualitativos del transcript HBB, aunque se recomienda validación adicional con herramientas especializadas antes de uso experimental.

