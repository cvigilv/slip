# Schüllerlab LIgand Priorization pipeline ("SLiP") - versión 1.0

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


