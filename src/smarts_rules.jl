const SMARTS_RULES = Dict(
 "glaxo" => "Glaxo",
 "dundee" => "dundee",
 "bms" => "BMS",
 "pains" => "PAINS",
 "surechembl" => "SureChEMBL",
 "mlsmr" => "MLSMR",
 "inpharmatica" => "Inpharmatica",
 "lint" => "LINT",
 "lilly" => "EliLilly",
 )

# Brenk filters (structural alerts for drug-likeness, 2008)
# Walters/Bruns filters (academic “bad actor” SMARTS sets)
# ToxAlerts / OECD QSAR Toolbox alerts (toxicophores, mutagenicity)
# REOS (Rapid Elimination of Swill) (Accelrys/BIOVIA HTS filters)
# Schrödinger SureFilter (commercial SMARTS set)
# ZINC Clean filters (used to sanitize ZINC database entries)

 #choose one filter
 #choose multiple filters
 #choose all filters
 #return names of filters (or dont)
