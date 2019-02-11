# Template - use the following (example for GBS):
# sed -e 's/__SPECIES__/GBS/g' browser-species-table.sql
# will still have to edit for columns
CREATE TABLE __SPECIES___BROWSER AS SELECT
  Tips.TIP AS number,
  Tips.Run AS run,
  CAST(CONCAT("Tip_", Tips.TIP) AS CHAR) AS treename,
  if (Studies.PrivacyCode = "Public" or
     (Studies.PrivacyCode IS NULL and
       (Tips.Run REGEXP '^[DES]RR' or Tips.Run REGEXP '^GC[AF]')),
     "Y", "N") AS public,
  IFNULL(__SPECIES__.DisplayName, IFNULL(Studies.SampleName, Tips.Run)) AS display_name,
  __SPECIES__.IsolationDate as final_year,
  if (MLST.MLSTDatabase = UPPER("__SPECIES__"), MLST, NULL) AS final_mlst,
  __SPECIES__.Serotype as final_serotype,
  __SPECIES__.Country as final_country,
  __SPECIES__.Host as final_host,
  __SPECIES__.Disease as final_disease,
  Studies.PrivacyCode as final_study,
  __SPECIES__.Comments as comments,
  ResistanceList.Aminoglycosides as aminoglycosides,
  ResistanceList.Betalactamases as betalactamases,
  ResistanceList.Colistin as colistin,
  ResistanceList.Fosfomycin as fosfomycin,
  ResistanceList.Fluoroquinolones as fluoroquinolones,
  ResistanceList.Glycopeptides as glycopeptides,
  ResistanceList.MLS as mls,
  ResistanceList.Phenicols as phenicols,
  ResistanceList.Rifampicin as rifampicin,
  ResistanceList.Sulfonamides as sulfonamides,
  ResistanceList.Tetracyclines as tetracyclines,
  ResistanceList.Trimethoprim as trimethoprim
FROM Tips
  LEFT JOIN Studies ON Tips.Run = Studies.Run
  LEFT JOIN MLST ON Tips.TIP = MLST.TIP
  LEFT JOIN ResistanceList ON Tips.TIP = ResistanceList.TIP
  RIGHT JOIN __SPECIES__ ON Tips.Run = __SPECIES__.Run
WHERE Tips.TIP IS NOT NULL;
UPDATE __SPECIES___BROWSER
  SET final_study = "Public"
  WHERE final_study IS NULL AND
        public = "Y";
UPDATE __SPECIES___BROWSER
  SET final_study = "slchen"
  WHERE final_study IS NULL;

