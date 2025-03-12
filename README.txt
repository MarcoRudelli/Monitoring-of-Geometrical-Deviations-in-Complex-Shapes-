
CONTENUTI
- Chipped1.obj’ e ‘Acceptable.obj’ (mesh modello IC e modello OC)
- ‘LinearFEM.m’ (funzione Matlab per il calcolo del LBS)
- ‘matlabmesh’ (toolkit utile per il caricamento delle mesh in Matlab)
- ‘generazione_mesh.m’ (script per la generazione delle mesh simulate)
- ‘remeshing_CAD_mesh.R’ (script per l’aumento del numero di triangoli della mesh IC originale)
- ‘calcolo_nearest_vertices.cpp’ (script per il calcolo dei vertici vicini ad ogni triangolo della mesh IC originale)
- ‘computazione_grafo.R’ (script per il calcolo del grafo utile al calcolo delle distanze geodetiche e per il calcolo dell’intorno di ogni triangolo)
- ‘calcolo_metrica_forward.cpp’ e ‘calcolo_metrica_backward_area.cpp’ (script per il caloclo delle metriche Forward, Backward e Area)
- ‘genera_statistiche_per_dfewma.R’ (script per lo smoothing geodetico delle metriche e la generazione dei csv utili alla sorveglianza delle metriche)
- ‘dfewma_personalizzato.cpp’ e ‘dfewma_classico.cpp’ (script per l’implementazione ad-hoc e classica della carta DFEWMA)
- ‘calcolo_risultati.R’ (script per mostrare i risultati della simulazione)



SETUP MATLAB:

Setup del toolkit ‘matlabmesh’:
	-> Set Path -> Add with Subfolders -> selezionare cartella ‘matlabmesh’ 

Algoritmo ICP:
	Add on:  https://www.mathworks.com/matlabcentral/fileexchange/12627-iterative-closest-point-method?s_tid=prof_contriblnk
	Nota: necessario aggiungere ‘res’ e ‘vi’ come ultimi due output aggiuntivi della funzione ‘icp’



ORDINE ESECUZIONE SCRIPT:

Eseguire ‘generazione_mesh.m’ per generare le mesh simulate e i risultati del LBS
Eseguire ‘remeshing_CAD_mesh.R’ per generare i file CSV per la mesh CAD con triangolazione originale e aumentata
Eseguire ‘calcolo_nearest_vertices.cpp’ per generare le liste dei vertici vicini ad ogni triangolo della triangolazione originale (utili per la definizione e l’utilizzo dei grafi per il calcolo delle distanze geodetiche) 
Eseguire ‘computazione_grafo.R’ per calcolare il grafo geodetico e l’insieme dei triangoli geodeticamente vicini ad ogni triangolo
Eseguire ‘calcolo_metrica_forward.cpp’ e ‘calcolo_metrica_backward_area.cpp’ per il calcolo delle metriche sulla superficie
Eseguire ‘genera_statistiche_per_dfewma.R’ per fare smoothing geodetico e generare i file per sorveglianza metriche
Eseguire ‘dfewma_personalizzato.cpp’ per sorveglianza DFEWMA metriche
Eseguire ‘dfewma_classico.cpp’ per sorveglianza DFEWMA LBS
Eseguire ‘calcolo_risultati.R’ per mostrare i risultati della simulazione
