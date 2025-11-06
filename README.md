
# Tarea 1: Análisis Funcional 


Entrega un script en Python que implemente de forma detallada un análisis funcional de los genes **COX4I2**, **ND1** y **ATP6**. El análisis debe utilizar librerías de Python para evaluar los procesos biológicos asociados a estos genes. Asegúrate de que tu código esté bien documentado y describa claramente los métodos y bases de datos utilizadas para obtener información funcional.

## Estructura del repositorio

```
/analisis-funcional/
├── data/
│   └── genes_input.txt        # Genes proporcionados para el análisis
├── scripts/
│   └── tu_script.py           # Tu script de análisis funcional
├── results/                   # Resultados generados por el script (opcional)
├── README.md                  # Este archivo
└── requirements.txt           # Dependencias del proyecto
```

## Instrucciones de entrega

- Haz un **fork** de este repositorio en tu cuenta de GitHub.
- Trabaja en tu fork y sube tu script a la carpeta `scripts/`.
- Tu script debe poder ejecutarse desde la línea de comandos (CLI), aceptando como mínimo el archivo de entrada y el archivo de salida como argumentos.
- No modifiques el archivo `genes_input.txt`.
- Documenta tu código explicando los métodos, librerías y bases de datos utilizadas.
- Puedes generar un archivo de resultados en la carpeta `results/` si lo consideras útil.

## Rúbrica de evaluación

La tarea se evaluará sobre un máximo de **10 puntos**, distribuidos según los siguientes criterios:

| Criterio | Descripción | Puntos |
|---------|-------------|--------|
| **1. Funcionalidad** | El script realiza correctamente el análisis funcional solicitado. | 4 |
| **2. Documentación** | El código está comentado y explica claramente los métodos y bases de datos utilizadas. | 2 |
| **3. Uso de librerías** | Se emplean librerías adecuadas para el análisis funcional (e.g., gprofiler-official, GOATOOLS, Enrichr, etc.). | 2 |
| **4. Formato y estilo** | El código sigue buenas prácticas de estilo y es legible. | 1 |
| **5. Automatización (CLI)** | El script acepta argumentos desde la línea de comandos. | 1 |

## Dependencias recomendadas

Incluye en `requirements.txt` las librerías necesarias para ejecutar tu script. Por ejemplo:

```
gprofiler-official
pandas
```

## Ejecución de la solución

* Instalación de dependencias: 
```
pip install -r requirements.txt
```
---

* Ejemplo de ejecución del script:
```
python .\scripts\script.py -i .\data\genes_input.txt --email usuario@correo.es -o resultados
```
---

* Parámetros de ejecución:

| Parámetro | Descripción | Valor por defecto | Ejemplo |
| --------------- | ---------------------------------------------------------------------------------- | ----------------- | ----------------------------- |
| `-i`, `--input` |  Ruta del archivo con los genes a analizar (en una sola línea, separados por comas) | — | `-i data/genes_input.txt` |
| `--email` | Correo electrónico obligatorio para las consultas a NCBI (Entrez) | — | `--email usuario@correo.es` |
| `-o`, `--out` | Nombre base para los archivos de salida (se añadirán sufijos de cada función) | `results` | `-o results` |

---

* Salida:
El script genera automáticamente la carpeta results/ (si no existe) y guarda los resultados de cada módulo con nombres coherentes:

| Archivo                      | Descripción                                                                  |
| ---------------------------- | ---------------------------------------------------------------------------- |
| `results/{out}_bio_python.txt`    | Anotaciones de los genes obtenidas con **Biopython (Entrez)**                |
| `results/{out}_go_enrichment.csv` | Resultados estructurados del análisis de enriquecimiento GO con **GOATOOLS** |
| `results/{out}_go_enrichment.txt` | Resumen legible del enriquecimiento GO                                       |
| `results/{out}_stringdb.csv`      | Resultados del análisis de enriquecimiento funcional con **STRINGdb**        |

---

