# Report di Splicing Alternativo da Annotazioni GTF

Il progetto analizza eventi di splicing alternativo nei geni. 
A partire da un file GTF filtrato su un singolo cromosoma, il notebook identifica e riassume eventi di splicing confrontando trascritti alternativi con uno di riferimento. 
I risultati sono esportati in un report strutturato.

## Descrizione
TEMA 3

Data l’annotazione GTF dei geni umani su un certo cromosoma (che si ottiene filtrando opportunamente l’annotazione completa Homo_sapiens.GRCh38.113.gtf.gz scaricabile da Ensemble Genome Browser),
produrre uno script/notebook che, per ogni gene annotato, selezioni un trascritto di riferimento e produca in output un report che riassuma gli eventi di splicing alternativo che gli altri trascritti 
del gene manifestano rispetto a quello preso come riferimento.

Precisamente, per ogni tipologia di evento, si deve riportare:

1) il numero di volte in cui si presenta globalmente nel cromosoma
2) la lista dei geni in cui si presenta
3) per ogni gene della lista del punto (2), la lista dei trascritti in cui l’evento si presenta

Le tipologie di evento da considerare sono exon skipping, 5’ competing, 3’ competing, intron retention e mutually exclusive exons. Si richiede di descrivere, per ogni tipologia, 
il criterio utilizzato per identificare l’evento.

NOTA BENE: non sono da prendere in considerazione eventi alternativi di inizio/fine della trascrizione.

AVVERTENZA: è necessario aggiungere il file *.gtf.gz nella cartella Annotazione_GTF.

## Caratteristiche
- Analisi di eventi di splicing alternativo: Exon Skipping, 5' Competing, 3' Competing, Intron Retention e Mutually Exclusive Exons.
- Creazione di report in Jupyter Notebook e un file Excel con dettagli sugli eventi di splicing.
- Supporta file GTF.

## Requisiti
- Python 3.x
- Librerie utilizzate:
  - pandas
  - tabulate
  - gffutils
  - bisect (standard)
  - csv (standard)
  - pathlib (standard)

## Struttura del Codice

Il progetto è composto da un unico script principale (`Progetto_Bio.ipynb`) che include:

### Funzioni principali

**save_gtf_and_create_db**(df_chr, output_dir, gtf_name, db_name) 
--> Salva il file GTF e crea un database utilizzando gffutils.

*Parametri*
df_chr (pandas.DataFrame) = DataFrame contenente le annotazioni GTF di un singolo cromosoma.
output_dir (Path) = Directory di destinazione in cui salvare il file .gtf e il database.
gtf_name (str) = Nome del file .gtf da salvare.
db_name (str) = Nome del file database .db da creare.

**get_reference_transcripts**(db) 
--> Estrae i trascritti di riferimento per un dato cromosoma dal database.

*Parametri*
db (gffutils.FeatureDB) = Database creato da un file GTF.

**get_exons**(db, transcript_id) 
--> Estrae gli esoni associati a un trascritto.

*Parametri*
db (gffutils.FeatureDB) = Database delle annotazioni GTF.
transcript_id (str) = Identificatore univoco di un trascritto.

**add_splicing_event**(event_type, gene_id, transcript_id) 
--> Aggiunge un evento di splicing alternativo alla lista degli eventi rilevati.

*Parametri*
event_type (str) = Tipo di evento di splicing.
gene_id (str) = ID del gene in cui si verifica l'evento.
transcript_id (str) = ID del trascritto che mostra l'evento.

**is_exon_skipped**(ref_exon, alt_exons) 
--> Verifica se un esone è stato "saltato" in un evento di exon skipping.

*Parametri*
ref_exon (gffutils.Feature) = Esone del trascritto di riferimento.
alt_exons (list of gffutils.Feature) = Lista di esoni del trascritto alternativo.

**overlap**(exon1, exon2) 
--> Calcola se due intervalli si sovrappongono.

*Parametri*
exon1 (gffutils.Feature) = Primo esone da confrontare.
exon2 (gffutils.Feature) = Secondo esone da confrontare.

**detect_exon_skipping**(gene_id, ref_middle_exons, trans_middle_exons, transcript_id) 
--> Rileva eventi di exon skipping.

*Parametri*
gene_id (str) = ID del gene analizzato.
ref_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto di riferimento.
trans_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto alternativo.
transcript_id (str) = ID del trascritto alternativo analizzato.

**detect_intron_retention**(gene_id, ref_exons, transcript_exons, transcript_id) 
--> Rileva eventi di intron retention.

*Parametri*
gene_id (str) = ID del gene analizzato.
ref_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto di riferimento.
trans_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto alternativo.
transcript_id (str) = ID del trascritto alternativo analizzato.

**detect_competing_sites**(gene_id, ref_middle_exons, trans_middle_exons, transcript_id) 
--> Rileva eventi di 5' e 3' competing.

*Parametri*
gene_id (str) = ID del gene analizzato.
ref_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto di riferimento.
trans_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto alternativo.
transcript_id (str) = ID del trascritto alternativo analizzato.

**detect_mutually_exclusive_exons**(gene_id, ref_middle_exons, trans_middle_exons, transcript_id) 
--> Rileva esoni mutuamente esclusivi.

*Parametri*
gene_id (str) = ID del gene analizzato.
ref_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto di riferimento.
trans_middle_exons (list of gffutils.Feature) = Esoni intermedi del trascritto alternativo.
transcript_id (str) = ID del trascritto alternativo analizzato.

### Strutture dati usate

**df** | pandas.DataFrame | DataFrame contenente l'annotazione GTF completa.

**df_cromosoma** | pandas.DataFrame | DataFrame filtrato per un cromosoma specifico.

**transcripts** | list → poi set | Lista (poi set) di trascritti estratti dal database.

**ref_transcripts** | set | Set dei trascritti di riferimento.

**trans_middle_exons_set** | set | Set degli esoni intermedi nei trascritti.

**transcript_exons_sorted** | list | Lista degli esoni di un trascritto ordinati.

**event_list** | list | Lista di eventi di splicing alternativo rilevati.

**db** | gffutils.FeatureDB | Database SQLite delle annotazioni GTF creato con gffutils.

**project_dir, output_dir** | pathlib.Path | Percorso al progetto e cartelle di output.

**splicing_event_counts** | dict (str → int) | Dizionario che tiene traccia del numero di eventi di splicing alternativo rilevati, suddivisi per tipo di evento.

**genes_per_event** | dict (str → set) | Dizionario che, per ogni tipo di evento di splicing, mantiene l'insieme (set) degli ID dei geni in cui quell'evento è stato osservato.

**transcripts_per_event** | dict (str → dict) | Dizionario che, per ogni tipo di evento di splicing, associa ad ogni gene (ID) i corrispondenti trascritti in cui l'evento è stato rilevato.
Praticamente: Tipo → Gene → Trascritti.

**exons_cache** | dict (str → list di esoni) | Memorizza in cache (per efficienza) gli esoni associati a ciascun trascritto, evitando di ricalcolarli più volte.
Gli esoni sono rappresentati come gffutils.Feature.

**transcripts_by_gene** | dict (str → list di str) | Raggruppa tutti i trascritti alternativi associati a ciascun gene (gene_id → lista di transcript_id).
