# Optimal Control of Acute Myeloid Leukemia (AML)
Questo repository contiene il progetto sviluppato per il corso di **Cibernetica Fisiologica** presso l'Università di Pisa (A.A. 2025/2026). 

L'obiettivo dello studio è l'applicazione della teoria del controllo ottimo per definire strategie di somministrazione chemioterapica che minimizzino il carico tumorale della Leucemia Mieloide Acuta, riducendo al contempo la tossicità per il paziente.

## Il Modello Matematico
La dinamica cellulare nel midollo osseo è stata simulata in **MATLAB** utilizzando un modello multi-compartimentale **SADLT** esteso, che modella:
* **S** (Staminali sane), **A** (Progenitrici), **D** (Differenziate)
* **L** (Staminali Leucemiche), **T** (Leucemiche Differenziate)
* **Competizione spaziale/nutrizionale** tra cellule sane e malate.
* **Risposta Immunitaria** modellata secondo la cinetica di Michaelis-Menten.

## Strategie di Controllo Implementate
Sono stati sviluppati e confrontati due approcci algoritmici per l'erogazione del farmaco chemioterapico:

### 1. Principio del Massimo di Pontryagin (PMP) - Open Loop
* **Controllo Continuo (Costo Quadratico):** Somministrazione fluida che bilancia eradicazione e tossicità.
* **Controllo Bang-Bang (Costo Lineare):** Alternanza tra dose massima e sospensione totale. Questa strategia fornisce una validazione matematica rigorosa alla pratica clinica delle *"Drug Holidays"* (vacanze terapeutiche).

### 2. Model Predictive Control (MPC) - Closed Loop
Poiché il PMP calcola la traiettoria ottima a priori (t=0), esso risulta "cieco" agli imprevisti clinici. Per garantire robustezza, è stato implementato un controllore MPC che misura lo stato del paziente a intervalli regolari e ricalcola l'ottimizzazione in tempo reale.
* **Risultato:** Simulando una ricaduta improvvisa (picco di cellule L), l'MPC rileva l'anomalia e adatta istantaneamente la dose, riuscendo a debellare la malattia dove il PMP fallisce.

##  Visualizzazione Risultati
<img width="1136" height="759" alt="confronto_controllo_continuo" src="https://github.com/user-attachments/assets/4141131d-6b6c-4be6-ad32-b9c8520a9a4c" />
<img width="1970" height="935" alt="drug_holidays" src="https://github.com/user-attachments/assets/8de22282-f9e0-466e-9647-33c66103899b" />
<img width="1035" height="938" alt="MPC_vs_PMP_con_disturbo" src="https://github.com/user-attachments/assets/963df686-1f51-4ff4-ad07-0d0e02bc85ee" />
