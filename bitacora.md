## 15 de Octubre de 2018
- Decido quitar la variable refAlign porque el archivo de entrada siempre tendrá el mismo nombre
- Agrego las variables de MoreFRags y posBInst
- Hago un commit en github en caso de que se dañe el inst2Bin dado que voy a quitar las cosas que no
  debo usar como Edis, y los archivos de salida que usaba Anibal como prueba.
## 1 de Noviembre de 2018
  - Tratar de poner a funcionar una versión secuencial
  [x] Hacer una copia de mappos, pasar dicha copia al radixsort y luego liberarla
  - Hacer un punto de control cuando se tenga una versión secuencial
  [x] En el archivo de alineamiento especial para el encoder debo agregar la suma total de errores que generó el ARG
  - Organizar el repositorio del ARG
  - Error Grave entre bases iguales Base1 0 Base2 0: Le debo prestar atención cuando la operacion sea s o S o i.
  - Revisar nuevamente la lectura de los datos desde el archivo de alineamiento
## 4 de Noviembre de 2018
  - El error que tenía era que estaba inicializando el puntero dentro del ciclo, entonces se reescribía con cada pasada del ciclo
  [] No olvidar liberar la memoria
  [] Aplicar lo del preámbulo
## 6 de Noviembre de 2018
  ### Objetivo: Encontrar porque el BinInst tiene como salida todos sus valores en 0
  - Del preámbulo están saliendo valores
  - Los elementos hasta totalReads tiene el valor del preámbulo