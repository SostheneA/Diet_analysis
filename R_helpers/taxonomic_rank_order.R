#
# https://api.gbif.org/v1/enumeration/basic/Rank
# https://en.wikipedia.org/wiki/Taxonomic_rank
#

ranknfile <- c("superdomain", "domain", "subdomain",
               "hyperkingdom", "superkingdom","kingdom","subkingdom","infrakingdom", "parvkingdom",
               "superphylum","phylum", "subphylum","infraphylum", "parvphylum",
               "gigaclass", "megaclass", "superclass","class","subclass","infraclass", "subterclass", "parvclass",
               "superdivision", "division", "subdivision", "infradivision",
               "superlegion","legion","sublegion","infralegion",
               "megacohort", "supercohort","cohort","subcohort","infracohort",
               "gigaorder", "magnorder","grandorder", "superorder","order","suborder","infraorder","parvorder",
               "section","subsection",
               "gigafamily", "superfamily","family","subfamily","infrafamily",
               "supertribe","tribe","subtribe","infratribe",
               "suprageneric_name","genus","subgenus","infragenus",
               
               "species_aggregate","species","infraspecific_name",
               "grex","subspecies","cultivar_group","convariety",
               "infrasubspecific_name","proles","race","natio","aberration","morph",
               "variety","subvariety","form","subform","pathovar","biovar","chemovar",
               "morphovar","phagovar","serovar","chemoform","forma_specialis",
               "cultivar","strain","other","unranked")

#ranknfile <- rev(ranknfile)