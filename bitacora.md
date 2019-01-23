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
## 10 de Noviembre
  ### Objetivo: Encontrar porque el BinInst tiene como salida todos sus valores en 0
  - Por cuestión de compatibilidad se trabaja a partir de la posición siguiente a la 0, por notación de alineadores
  #### Definición de los campos:
    8 bits que contienen un PRELUDE (hasta el momento, se va a cambiar por que en 8 bits quepan los preludios de dos reads (4 bits por preludio))
    8 bits los menos significativos del offset
    8 bits que se dividen así:
      2 bits que faltan del offset
      1 bits que indica si este es el último error
      3 bits opcode (codigo de la operación)
      2 bits que indican la diferencia entre las bases

  EUREKA
## 19 de Noviembre
  ### Objetivo: Definir el nuevo preludio
  - Definición: Si mal no recuerdo, la idea es sacar los preludios del BinInst y ponerlos en su propio
  arreglo de preludios, en general, por cada 8 bits caben 2 preludios
  - Defino un arreglo de preludios:
    - Recordar que un preludio es un conjunto de 4 bits, el primer bit indica si el siguiente read mapea en la misma posición, los siguientes tres bits corresponden a el código del matching.
    - El arreglo de preludios debe ser de tamaño TotalRead/2, en caso de que el total de reads sea impar,
    entonces se define dicho tamaño como floor(TotalReads/2)+1

### 23 de Noviembre
  ### Objetivo: Comenzar con el paralelismo
  - Planificación ( Distribución de la carga, balanceo de carga dinámica y balanceo de carga estática )
  - En OpenMp está el parámetro schedule
  _PROXIMAS ACTIVIDADES_
  1. HACER PRUEBAS DE QUE LAS COSAS CORRAN CORRECTAMENTE (SEMÁNTICA)

  2. Tiempo esperado, de 2 a tres veces mejor 
  3. Mirar el tema de la distribucion de la carga de trabajo en los reads usando los dos modelos
  4. Verificar los preludios
  
  (Modelo estático: Todos los hilos tienen el mismo número de reads excepto el ultimo)
  (Modelo dinámico: los hilos pelean por chunks, tema de las condiciones de carrera )

### 29 de Noviembre
  ### Objetivo: Extraer y probar la versión secuencial
  - Antes de comenzar la paralelización, organizo la versión secuencial


### RECESO
  ### MOTIVO: Daño del computador pospuso la actividad

# 2019

### 10 de Enero
  ### Objetivo: Retomar el trabajo de grado
  #### Estrategia:
    -[x] Organizar el repositorio de github y acomodarlo a la estrategia de desarrollo de git-flow
    -[] Finalizar de una vez por todas la versión secuencial, garantizando que funciona para las diferentes cantidades de datos posibles.
    -[] Comenzar a organizar la documentación.
  ### Puntos clave:
    - El nombre del archivo que se va a procesar esta quemado en el fopen, se debe pasar este parámetro por la línea de comandos.
    - Parece que el ARG no generó los 200.000 reads.
    - Descargue el ARG_ENCODER y volví a generar un archivo de alineamiento, verifiqué que tuviera los 200000 reads, y lo paśe al encoder, en esta ocasión funcionó correctamente. De aquí surge un problema, es probable que bajo ciertas circunstancias la lectura del archivo de alineamiento fracase, hay que verificar dichas circunstancias.

### 17 de Enero
  ### Puntos clave:
    - Se hacen pruebas de los preambulos en la versión secuencial

### 22 de Enero
  ### Puntos Clave:
  ![ESTRUCTURAS](/IMG/ESTRUCTURAS.jpeg)
  - 

