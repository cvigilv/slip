# Formato archivos requeridos

El formato para los 3 archivos necesitados para correr _SLiP_ se encuentra a continuación. Cada entrada se organiza de la siguiente manera:

- El nombre asignado para ese archivo en el archivo de configuraciones.

- Una lista enumerada que corresponde a el tipo de informacion que debe contener dicha columna en el archivo.

Se debe tomar en cuenta que el separador de estos archivos debe ser *tab* ('\t'), de lo contrario arrojará error.

#### _"Interactions file"_ o salida de _pipeline_ para predicir interacciones proteína-ligando.

1. _Fold_
2. _"Query ligand ChEMBL ID'_ / Identificador ChEMBL para el ligando consulta.
3. _"Hit target ChEMBL ID"_ / Identificador ChEMBL para el blanco predicho.
4. Nombre de medida de similitud (ejemplo _Tc_, _LR_, _Srel_, etc.)) 
5. _"Hit ligand ID"_ / Identificador para ligando encontrado.
6. _"Query target ID"_ / Identificador para blanco consulta.
7. _TP_

#### _"Broad file"_ o selección amplia de la base de datos ChEMBL.

1. Identificador ChEMBL para el blanco.
2. Identificador ChEMBL para el ligando.
5. Identificador(es) Pfam para el blanco.
-3. Fase clínica máxima para el ligando.
-2. Número de átomos pesados del ligando.


#### _"SMILES file"_ o archivo relacional "SMILES canonicos - ChEMBL ID ligandos".

1. SMILES canónicos.
2. ChEMBL ID ligando.

