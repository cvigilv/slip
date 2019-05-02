# Schüllerlab LIgand Priorization pipeline ("SLiP") - versión 1.0

Tabla de contenido
-----
* [¿Qué es SLiP?](#qué-es-slip)
* [Dependencias](#dependencias)
* [¿Cómo utilizarlo?](#cómo-utilizarlo)
* [Configuraciones](#configuraciones)
	* [Sección <em>[Input]</em>](#sección-input)
		* [<em>Interactions file:</em>](#interactions-file)
		* [<em>Broad file</em>:](#broad-file)
  		* [<em>ChEMBL version</em>:](#chembl-version)
 	* [Sección <em>[Options]</em>](#sección-options)
  		* [<em>Keep target</em>:](#keep-target)
  		* [<em>Keep ligand</em>:](#keep-ligand)
  		* [<em>Pfam cutoff</em>:](#pfam-cutoff)
  		* [<em>Top X entries</em>:](#top-x-entries)
  		* [<em>Similarity measure used</em>:](#similarity-measure-used)
  		* [<em>min(Similarity)</em>:](#minsimilarity)
  		* [<em>max(Similarity)</em>:](#maxsimilarity)
  		* [<em>min(Clinical Phase)</em>:](#minclinical-phase)
  		* [<em>max(Clinical Phase)</em>:](#maxclinical-phase)
 	* [Sección <em>[Output]</em>](#sección-output)
	 	* [<em>Output directory</em>:](#output-directory)
  		* [<em>Output file</em>:](#output-file)
  		* [<em>Generate plots</em>:](#generate-plots)
  		* [<em>Prepare file</em>:](#prepare-file)


## ¿Qué es SLiP?

SLiP o _"Schüllerlab Ligand Priorization pipeline"_ corresponde a una herramienta escrita en Python que tiene por finalidad filtrar y priorizar el archivo de salida proveniente del _pipeline_ para predecir interacciones proteina-ligando; la cual puede emplearse en análisis posteriores, como simulaciones de _docking_ molecular o visualización de predicciones interesantes.

## Dependencias

SLiP necesita de los siguientes archivos para poder correr:

- Archivo de salida del _pipeline_ de predicciones proteína-ligando.

- Archivo con una selección amplia o _"broad"_ de ChEMBL.

- Archivo con SMILES para los ligandos del archivo de salida.

- Archivo de configuraciones para el programa.

Para conocer el formato de estos archivos, referir al archivo [FORMAT.md](FORMAT.md) dentro del repositorio.

Además SliP necesita:

- Python version >= 3.0.

- pandas

- matplotlib

- seaborn

- pymysql

## ¿Cómo utilizarlo?

En la linea de comando, se puede correr SLiP de la siguiente manera:

        $ python3 slip.py slip.config

## Configuraciones

_SLiP_ cuenta con un [archivo de configuraciones](slip.ini) que le permite modificar ciertas configuraciones estándar con las que corre el programa.

### Sección _[Input]_

#### _Interactions file:_ 
Corresponde al archivo de salida del pipeline de predicciones proteína-ligando. El formato permitido debe ser el siguiente:

| Número columna | Contenido columna | Python type() | Comentario |
|:--------------:|:------------------|:-------------:|:-----------|
| 1 | _Fold_| int | |
| 2 | _Query Ligand ChEMBL ID_ | str | |
| 3 | _Hit Target ChEMBL ID_ | str | |
| 4 | _Similarity measure_ | float | Debe estar en un rango de 0.0 a 1.0 (normalizado) |
| 5 | _Query Target ID_ | str | Columna ignorada |
| 6 | _Hit Ligand ID_  | str | **Debe ser ChEMBL ID** |
| 7 | _TP_ | int | Debe ser 0 o 1 |

#### _Broad file_: 
Archivo de interacciones proteína-ligando amplia, referirse a [FORMAT.md](FORMAT.md) para más información.

#### _ChEMBL version_: 
Versión de ChEMBL utilizada para la validación temporal y contra otras especies para eliminar casos que se han probado en otras especies o que se encuentran en versiones de ChEMBL más acutalizadas.

### Sección _[Options]_

#### _Keep target_:
Esta opción permite generar un archivo de salida que solo contanga predicciones para un blanco especificado. Esta entrada debe ser un ChEMBL ID para el blanco que se desea analizar.

#### _Keep ligand_: 
Esta opción permite generar un archivo de salida que solo contenga predicciones para un ligando especifico. Esta entrada debe ser un ChEMBL ID para el ligando que se desea analizar.

#### _Pfam cutoff_: 
Número de Pfam ID's que pueden compartir el ligando consultado y el blanco predicho. De ser igual o mayor a este valor, se eliminará la entrada.

#### _Top X entries_: 
Dejar el top X de entradas según el valor del _similarity measure_.

#### _Similarity measure used_: 
Medida de similitud o de puntuación de la predicción. Se utiliza para generar gráficos y tablas acorde a esta medida.

#### _min(Similarity)_: 
Limite inferior inclusivo utilizado para filtrar por similitud. Debe ser un valor entre 0 y 1.

#### _max(Similarity)_: 
Limite superior no inclusivo utilizado para filtrar por similitud. Debe ser un valor entre 0 y 1.

#### _min(Clinical Phase)_: 
Limite inferior inclusivo utilizado para filtrar por fase clínica del ligando. Debe ser un valor entre 0 y 4.

#### _max(Clinical Phase)_: 
Limite superior no inclusivo utilizado para filtrar por la fase clínica del ligando. Debe ser un valor entre 0 y 4.

### Sección _[Output]_

#### _Output directory_: 
Directorio donde se almacenaran todos los archivos de salida de SLiP. De dejarse vacío, se creará una carpeta con el nombre del archivo de interacciones con la terminación \*-slip.

#### _Output file_: 
Nombre del archivo de salida filtrado. De dejarse vacío, se colocará la terminación \*-slip al archivo.

#### _Generate plots_: 
Genera gráficos de distribuciones a partir de cada paso de filtrado. Debe ser "True" o "False".

#### _Prepare file_: 
Preparar el archivo de salida utilizando MOE y el script [slip\_prepare.svl](slip\_prepare.svl). Debe ser "True" o "False".
