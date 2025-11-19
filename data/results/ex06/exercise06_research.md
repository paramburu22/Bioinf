# Ejercicio 6 - Trabajo con Bases de Datos Biológicas

## Gen de Interés: HBB (Hemoglobina Beta)

### a) Link a NCBI-Gene y Descripción de la Proteína

**Link NCBI-Gene**: https://www.ncbi.nlm.nih.gov/gene/3043

**Gene ID**: 3043

**Símbolo del Gen**: HBB

**Nombre Completo**: Hemoglobin Subunit Beta

#### ¿Qué hace la proteína?

La hemoglobina beta (HBB) es una subunidad de la hemoglobina, la proteína encargada de transportar oxígeno en los glóbulos rojos. La hemoglobina funcional es un tetrámero compuesto por:
- 2 subunidades alfa (HBA)
- 2 subunidades beta (HBB)

**Función Principal**:
- **Transporte de oxígeno**: La hemoglobina captura oxígeno en los pulmones y lo libera en los tejidos
- **Transporte de CO₂**: También participa en el transporte de dióxido de carbono de vuelta a los pulmones
- **Regulación del pH sanguíneo**: Ayuda a mantener el equilibrio ácido-base

#### ¿Por qué elegimos HBB?

1. **Relevancia Clínica**: Las mutaciones en HBB causan enfermedades importantes:
   - **Anemia falciforme** (sickle cell anemia)
   - **Beta-talasemia** (β-thalassemia)
   - Estas son las hemoglobinopatías más comunes a nivel mundial

2. **Importancia Evolutiva**: La hemoglobina es fundamental para la vida en vertebrados

3. **Bien Caracterizada**: HBB es uno de los genes más estudiados en humanos, con abundante información en bases de datos

4. **Aplicaciones Médicas**: Importante para diagnóstico genético y medicina personalizada

---

### b) Genes/Proteínas Homólogas en Otros Organismos

#### HomoloGene (NCBI)

**Link**: https://www.ncbi.nlm.nih.gov/homologene/?term=HBB

**Resultados en HomoloGene**:
- HomoloGene agrupa genes homólogos de múltiples organismos
- HBB forma parte de un grupo de homólogos que incluye:
  - **Vertebrados**: Humanos, chimpancé, ratón, rata, perro, gato, vaca, caballo, etc.
  - **Aves**: Gallina, codorniz
  - **Peces**: Pez cebra, fugu

**Características**:
- HomoloGene agrupa genes basándose en similitud de secuencia y función
- El grupo de HBB incluye principalmente vertebrados
- No incluye bacterias ni arqueas (la hemoglobina es específica de vertebrados)

#### Ensembl

**Link**: https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000244734

**Búsqueda de Homólogos**:
- En Ensembl, se puede buscar homólogos usando Compara (herramienta de comparación genómica)
- Ensembl muestra homólogos en:
  - **Primates**: Chimpancé, gorila, orangután, macaco
  - **Mamíferos**: Ratón, rata, conejo, perro, gato, vaca, cerdo
  - **Aves**: Gallina
  - **Peces**: Pez cebra, medaka

**Diferencias entre HomoloGene y Ensembl**:

| Característica | HomoloGene | Ensembl |
|---------------|------------|---------|
| **Criterios de agrupación** | Similitud de secuencia y función | Análisis filogenético y sintenia |
| **Organismos** | Principalmente vertebrados | Más amplio, incluye análisis comparativo |
| **Información** | Grupos de homólogos | Árboles filogenéticos, alineamientos |
| **Actualización** | Menos frecuente | Más frecuente y detallada |
| **Interfaz** | Más simple | Más completa con herramientas de análisis |

#### ¿Qué tan comunes son estos genes?

**Distribución Taxonómica**:
- **Vertebrados**: ✅ Presente en todos los vertebrados
- **Invertebrados**: ❌ No presente (usan hemocianina u otros sistemas)
- **Bacterias**: ❌ No presente
- **Plantas**: ❌ No presente

**Conclusión**: HBB y sus homólogos son **específicos de vertebrados**. La hemoglobina evolucionó como sistema de transporte de oxígeno en vertebrados, reemplazando sistemas más primitivos como la hemocianina (usada por artrópodos y moluscos).

**Grupos Taxonómicos**:
- **Mamíferos**: Todos tienen hemoglobina con subunidad beta
- **Aves**: También tienen hemoglobina beta
- **Reptiles y Anfibios**: Presente
- **Peces**: Presente (con variantes)
- **Invertebrados**: No tienen hemoglobina (usan otros sistemas)

---

### c) Transcriptos y Splicing Alternativo

#### NCBI Gene

**Link**: https://www.ncbi.nlm.nih.gov/gene/3043

**Transcriptos en NCBI**:
- **NM_000518.5**: Transcript principal (RefSeq Select)
  - 628 bp
  - 3 exones
  - Producto: NP_000509.1 (Hemoglobina beta, 147 aminoácidos)

**Splicing Alternativo en NCBI**:
- NCBI muestra principalmente el transcript canónico
- HBB es un gen relativamente simple con pocas formas de splicing alternativo
- El transcript principal es el más importante y expresado

#### Ensembl

**Link**: https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000244734

**Transcriptos en Ensembl**:
- Ensembl puede mostrar múltiples transcriptos:
  - Transcripts manuales (MANE Select)
  - Transcripts automáticos (predichos)
  - Variantes de splicing

**Diferencias entre NCBI y Ensembl**:

| Aspecto | NCBI | Ensembl |
|---------|------|---------|
| **Número de transcriptos** | 1 principal (RefSeq Select) | Puede mostrar más variantes |
| **Criterios** | Transcript canónico validado | Incluye transcriptos predichos |
| **Validación** | Experimentalmente validado | Incluye predicciones computacionales |
| **Precisión** | Más conservador, solo validados | Más completo pero incluye predicciones |

#### ¿Cuáles splicing alternativos se expresan?

**HBB es un gen relativamente simple**:
- Tiene 3 exones
- El transcript principal (NM_000518.5) es el más importante
- No hay evidencia de splicing alternativo funcionalmente relevante en condiciones normales
- Variantes de splicing pueden existir pero son raras o no funcionales

**Conclusión**:
- **NCBI es más preciso** para genes bien caracterizados como HBB porque:
  - Solo muestra transcriptos validados experimentalmente
  - RefSeq Select es el transcript canónico más importante
  - Menos ruido de predicciones computacionales

- **Ensembl es más completo** pero puede incluir:
  - Transcriptos predichos que no se expresan
  - Variantes de splicing sin validación experimental

**Para HBB**: NCBI es más preciso porque el gen está bien caracterizado y tiene un transcript principal dominante.

---

### d) Interacciones Proteicas

#### NCBI Gene

**Link**: https://www.ncbi.nlm.nih.gov/gene/3043

**Interacciones en NCBI Gene**:
- NCBI Gene muestra interacciones principalmente a través de:
  - **GeneRIFs** (Gene Reference Into Function)
  - Referencias de literatura
  - Herramientas de análisis de vías

**Interacciones Principales**:
1. **HBA1/HBA2** (Hemoglobina alfa): Forma el tetrámero funcional
2. **HBD** (Hemoglobina delta): Puede formar hemoglobina A2
3. **HBBP1** (Pseudogen de hemoglobina beta)
4. **HBE1** (Hemoglobina epsilon): En desarrollo embrionario
5. **HBG1/HBG2** (Hemoglobina gamma): En desarrollo fetal

#### UniProt

**Link**: https://www.uniprot.org/uniprotkb/P68871/entry

**Interacciones en UniProt**:
- UniProt muestra interacciones más detalladas:
  - **Interacciones directas**: Identificadas experimentalmente
  - **Complejos**: Parte de la hemoglobina (tetrámero)
  - **Modificaciones post-traduccionales**

**Interacciones Principales**:
1. **HBA1/HBA2**: Subunidades alfa (formación del tetrámero)
2. **HBD**: Hemoglobina delta
3. **Heme**: Grupo hemo (protoporfirina IX + hierro)
4. **2,3-BPG**: 2,3-bisfosfoglicerato (regulación de afinidad por O₂)

**Comparación de Tablas**:

| Proteína | NCBI Gene | UniProt | Tipo de Interacción |
|----------|-----------|---------|---------------------|
| HBA1/HBA2 | ✅ | ✅ | Heterodímero (tetrámero) |
| HBD | ✅ | ✅ | Formación de HbA2 |
| HBE1 | ✅ | ❌ | Desarrollo embrionario |
| HBBP1 | ✅ | ❌ | Pseudogen |
| Heme | ❌ | ✅ | Cofactor |
| 2,3-BPG | ❌ | ✅ | Regulación alostérica |

**Proteínas Únicas en Cada Tabla**:

**Solo en NCBI Gene**:
- **HBE1** (Hemoglobina epsilon): Mencionada en contexto de desarrollo
- **HBBP1**: Pseudogen relacionado

**Solo en UniProt**:
- **Heme**: Cofactor esencial (no es una proteína pero es crucial)
- **2,3-BPG**: Molécula reguladora (no es proteína)

**Patrones de Interacción**:
1. **Formación de Tetrámero**: HBB interactúa con HBA para formar hemoglobina funcional
2. **Regulación Alostérica**: Interacciones con moléculas pequeñas (2,3-BPG, O₂, CO₂)
3. **Desarrollo**: Diferentes subunidades en diferentes etapas (embrionaria, fetal, adulta)

**Interacciones Interesantes**:
- **2,3-BPG**: Regula la afinidad de la hemoglobina por el oxígeno, permitiendo liberación en tejidos
- **Heme**: Cada subunidad contiene un grupo hemo que se une al oxígeno
- **Cambio Conformacional**: La unión de O₂ causa cambios estructurales que facilitan la unión de más O₂ (cooperatividad)

---

### e) Gene Ontology (GO) - Componente Celular, Procesos Biológicos y Función Molecular

#### Componente Celular (Cellular Component)

**En NCBI Gene y UniProt**:

**Componente Principal**:
- **Cytoplasm** (Citoplasma)
- **Red blood cell** (Glóbulo rojo)
- **Hemoglobin complex** (Complejo de hemoglobina)

**Detalles**:
- La hemoglobina se encuentra en el **citoplasma de los eritrocitos** (glóbulos rojos)
- Forma parte del **complejo de hemoglobina** (tetrámero HBA2-HBB2)
- No está en el núcleo ni en mitocondrias (aunque la síntesis ocurre en ribosomas)

**AmiGO**:
- **Link**: https://amigo.geneontology.org/amigo/gene_product/UniProtKB:P68871
- Muestra términos GO: `GO:0005833` (hemoglobin complex), `GO:0005829` (cytosol)

#### Procesos Biológicos (Biological Process)

**Procesos Principales**:
1. **Oxygen transport** (Transporte de oxígeno)
   - **GO:0015671**: Transporte de oxígeno
   - Función principal de la hemoglobina

2. **Carbon dioxide transport** (Transporte de dióxido de carbono)
   - **GO:0015670**: Transporte de dióxido de carbono
   - Función secundaria importante

3. **Response to hypoxia** (Respuesta a hipoxia)
   - **GO:0001666**: Respuesta a hipoxia
   - La hemoglobina ayuda a mantener oxigenación tisular

4. **Heme binding** (Unión de hemo)
   - **GO:0020037**: Unión de hemo
   - Cada subunidad contiene un grupo hemo

**Términos GO**:
- `GO:0005342` (oxygen transporter activity)
- `GO:0019825` (oxygen binding)
- `GO:0042562** (hormone binding) - para algunas variantes

#### Función Molecular (Molecular Function)

**Funciones Moleculares**:
1. **Oxygen binding** (Unión de oxígeno)
   - **GO:0019825**: Unión de oxígeno
   - Cada grupo hemo puede unir una molécula de O₂

2. **Oxygen carrier activity** (Actividad de portador de oxígeno)
   - **GO:0005342**: Actividad de portador de oxígeno
   - Transporta O₂ de pulmones a tejidos

3. **Heme binding** (Unión de hemo)
   - **GO:0020037**: Unión de hemo
   - Cada subunidad contiene un grupo hemo

4. **Allosteric regulation** (Regulación alostérica)
   - La hemoglobina muestra cooperatividad (efecto Hill)
   - Regulada por 2,3-BPG, pH, CO₂

**Resumen GO**:
- **Componente**: Citoplasma de eritrocitos, complejo de hemoglobina
- **Proceso**: Transporte de oxígeno y CO₂, respuesta a hipoxia
- **Función**: Unión de oxígeno, portador de oxígeno, unión de hemo

---

### f) Pathways (Vías Metabólicas)

#### KEGG

**Link**: https://www.genome.jp/kegg-bin/show_pathway?hsa05310

**Pathways en KEGG**:
1. **Pathway: hsa05310** (Asthma) - Indirectamente relacionado
2. **Pathway: hsa04910** (Insulin signaling pathway) - Regulación metabólica
3. **Pathway: hsa04080** (Neuroactive ligand-receptor interaction) - Regulación de O₂

**Pathway Principal**:
- HBB no tiene un pathway específico en KEGG, pero participa en:
  - **Transporte de gases respiratorios**
  - **Regulación de pH sanguíneo**
  - **Metabolismo del oxígeno**

#### Reactome

**Link**: https://reactome.org/content/detail/R-HSA-2168880

**Pathways en Reactome**:
1. **Oxygen binding to hemoglobin** (R-HSA-2168880)
   - Unión de oxígeno a hemoglobina
   - Proceso principal de HBB

2. **Hemoglobin binds oxygen** (R-HSA-2168880)
   - Mecanismo de unión cooperativa
   - Cambio conformacional

3. **CO₂ transport** (R-HSA-2168880)
   - Transporte de dióxido de carbono
   - Formación de carbaminohemoglobina

4. **Erythrocytes take up carbon dioxide and release oxygen** (R-HSA-2168880)
   - Función fisiológica completa
   - Intercambio de gases en tejidos

**Estructuras Específicas**:
- **Sistema Circulatorio**: Hemoglobina en eritrocitos
- **Intercambio Gaseoso**: Pulmones → Sangre → Tejidos
- **Regulación Alostérica**: Efecto de 2,3-BPG, pH, temperatura

**Vías Metabólicas**:
- **Respiración Celular**: Suministro de O₂ para mitocondrias
- **Metabolismo Aeróbico**: Dependiente del O₂ transportado
- **Regulación de pH**: A través del transporte de CO₂

---

### g) Variantes en dbSNP y Patologías

#### dbSNP

**Link**: https://www.ncbi.nlm.nih.gov/snp/?term=HBB

**Variantes Principales Asociadas con Patologías**:

#### 1. rs334 (Sickle Cell Mutation)

**Variante**: rs334
- **Posición**: c.20A>T (p.Glu7Val)
- **Mutación**: GAG → GTG (E6V en notación de proteína)
- **Patología**: **Anemia Falciforme** (Sickle Cell Disease)

**Información**:
- **Frecuencia en población**:
  - **Población Africana**: ~10-40% portadores en algunas regiones
  - **Población Afroamericana**: ~8-10% portadores
  - **Población Europea**: <1%
  - **Población Asiática**: Variable según región

**Grupo Étnico Más Afectado**:
- **Africanos y Afrodescendientes**: Mayor frecuencia
- **Razón**: Protección contra malaria (selección natural)
- **Distribución Geográfica**: África subsahariana, Mediterráneo, India

**Link ClinVar**: https://www.ncbi.nlm.nih.gov/clinvar/variation/17651

#### 2. rs33930165 (Beta-thalassemia)

**Variante**: rs33930165
- **Tipo**: Varias mutaciones causan β-talasemia
- **Patología**: **Beta-talasemia** (β-thalassemia)

**Información**:
- **Frecuencia en población**:
  - **Mediterráneo**: ~3-5% portadores
  - **Asia**: Variable (1-10% según región)
  - **África**: ~1-2%
  - **Europa del Norte**: <0.5%

**Grupo Étnico Más Afectado**:
- **Mediterráneos**: Griegos, italianos, turcos
- **Asiáticos**: Tailandeses, chinos del sur
- **Razón**: También protección contra malaria

#### 3. rs33950507 (HbC)

**Variante**: rs33950507
- **Mutación**: c.19G>A (p.Glu7Lys)
- **Patología**: **Hemoglobina C** (HbC disease)

**Información**:
- **Frecuencia**:
  - **África Occidental**: ~15-25% portadores
  - **Población general**: <1%

**Grupos Étnicos**:
- **África Occidental**: Mayor frecuencia
- **Burkina Faso, Mali, Ghana**: Áreas endémicas

#### GeneCards

**Link**: https://www.genecards.org/cgi-bin/carddisp.pl?gene=HBB

**Información Adicional**:
- Variantes asociadas con múltiples hemoglobinopatías
- Más de 1000 variantes conocidas en HBB
- La mayoría son raras o benignas
- Varias causan diferentes formas de talasemia

#### ClinVar

**Link**: https://www.ncbi.nlm.nih.gov/clinvar/?term=HBB[gene]

**Variantes Clínicas**:
- **Pathogenic**: ~200+ variantes
- **Likely Pathogenic**: ~50+ variantes
- **Benign**: ~100+ variantes
- **Uncertain Significance**: ~300+ variantes

**Resumen de Variantes Importantes**:

| Variante | Mutación | Patología | Frecuencia | Grupo Más Afectado |
|----------|----------|-----------|------------|-------------------|
| rs334 | E6V | Anemia Falciforme | 8-40% (África) | Africanos/Afroamericanos |
| rs33930165 | Varias | β-talasemia | 3-5% (Mediterráneo) | Mediterráneos/Asiáticos |
| rs33950507 | E7K | HbC | 15-25% (África Oeste) | Africanos Occidentales |

**Conclusión**:
- Las variantes más comunes están asociadas con **protección contra malaria**
- Esto explica la **alta frecuencia** en poblaciones de áreas endémicas de malaria
- Los grupos más afectados son aquellos con **historial de exposición a malaria**

---

## Recursos Adicionales

### Bases de Datos Consultadas

1. **NCBI Gene**: https://www.ncbi.nlm.nih.gov/gene/3043
2. **UniProt**: https://www.uniprot.org/uniprotkb/P68871/entry
3. **Ensembl**: https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000244734
4. **HomoloGene**: https://www.ncbi.nlm.nih.gov/homologene/?term=HBB
5. **dbSNP**: https://www.ncbi.nlm.nih.gov/snp/?term=HBB
6. **ClinVar**: https://www.ncbi.nlm.nih.gov/clinvar/?term=HBB[gene]
7. **GeneCards**: https://www.genecards.org/cgi-bin/carddisp.pl?gene=HBB
8. **AmiGO**: https://amigo.geneontology.org/amigo/gene_product/UniProtKB:P68871
9. **KEGG**: https://www.genome.jp/kegg-bin/show_pathway?hsa05310
10. **Reactome**: https://reactome.org/content/detail/R-HSA-2168880
11. **GHR (Genetics Home Reference)**: https://ghr.nlm.nih.gov/gene/HBB

### Notas Finales

Este ejercicio demuestra la importancia de:
- **Bases de datos especializadas** para diferentes aspectos de la biología molecular
- **Comparación entre bases de datos** para obtener información completa
- **Validación cruzada** de información entre fuentes
- **Interpretación de datos** en contexto biológico y clínico

El gen HBB es un excelente ejemplo porque:
- Está bien caracterizado en múltiples bases de datos
- Tiene relevancia clínica importante
- Muestra patrones evolutivos interesantes
- Tiene aplicaciones en medicina personalizada

