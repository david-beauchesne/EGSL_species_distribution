# Run init.r before other scripts (set the working directory and load the pkges )
rm(list=ls())

# -----------------------------------------------------------------------------
# PROJECT:
#    Evaluating the structure of the communities of the estuary
#       and gulf of St.Lawrence
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# DETAILS:
#   The goal of this script is to run JSDMs from the HMSC package on data from
#       the annual northern Gulf DFO trawl survey
# -----------------------------------------------------------------------------



# Importing data required data from RawData
    northPluri <- readRDS('../../../PhD_RawData/data/NordGolfe_RelevePluriSp_MPO/northPluri.rds')
    northPluri[,'N_EspSci'] <- str_trim(northPluri[,'N_EspSci'], side = 'both')
    northPluri[,'N_EspF'] <- as.character(northPluri[,'N_EspF'])

# Species name list & number of records per year
    # Count number of observations per species
    countRec <- numeric(nrow(northPluri))
    spNames <- aggregate(countRec ~ EspGen + N_EspSci + N_EspF, data = cbind(northPluri, countRec), FUN = length) %>%
                    .[order(.[, 'EspGen']), ]

# Which taxa should be combined?
    # Form groups of taxa as a vector, with taxa to combine seperated with ' | '
    # Main taxa that will be used for aggregation should be place first

    # Taxa are kept, to the extent possible, at the species level. Groups that are too coarse are generally removed from the dataset. Combined taxa are generally species that are closely related with similar functional roles in ecosystems and that are hard to distinguish in the field. For example, species of the genus Buccinum are aggregated as they have similar functional roles and can be very hard to differentiate. Certain large groups like Porifera were also retained, as they represent a widely distributed group and their identification in the field is hard to accomplish. Removing Porifera in favour of species would likely result in a gross underestimates of this group in the St. Lawrence and likely include misidentifications.

    # Every removal or aggregation is documented in the following code to allow for increased transparency and reproducibility.

    # Additionnal notes:
        # 'Bryozoa | Alcyonidium sp. | Alcyonidium pachydermatum | Reteporella grimaldii | Securiflustra securifrons | Caberea ellisii' # Bryozoa retained, coarse but represent a widely distributed group that would be vastly underestimated if left out. Other 4 species and 1 genus are larger, errected species more easily recognized, hence they are kept individually. Since there is not a lot of data to describe them, however, it may be wise to group them all under the bryozoa umbrella for the analysis.

    # To verify in particular:
        # - Myctophidae
        # - Gaidropsarus sp.
        # - Artediellus sp.
        # - Liparidae (stays at the species level)
        # - Lycodes (stays at the species level) + Lycodes terraenovae & Lycodes esmarkii
        # - Lycenchelys (stays at the species level)
        # - Adjust Eumicrotremus spinosus

    combineTaxa <- c('Strongylocentrotus sp. | Strongylocentrotus droebachiensis', # Identification: hard to distinguish
                     'Sebastes sp. | Sebastes norvegicus', # Identification: hard to distinguish
                     'Porifera | Tentorium semisuberites | Radiella hemisphaerica | Polymastia sp.', # Combine, hard to distinguish in the field. Porifera is coarse, but represent a widely distributed group that would be vastly underestimated if left out. We therefore combine them at the 'Porifera' level, only leaving 'Stylocordyla borealis' as a species, as it is easily identified
                     'Atolla wyvillei | Atolla sp.', # Very likely the same species
                     'Polynoidae | Euphrosine borealis | Harmothoe sp. | Eunoe nodosa | Laetmonice filicornis', # Identification: can be hard to distinguish in the field and they are species that generally have similar functional roles and habitat requirements. Only 'Aphroditella hastata' is left as a species, as it is frequently caught and easily identifiable due to its size
                     'Ascidiacea | Pelonaia corrugata | Botrylloides sp. | Ascidia sp. | Cnemidocarpa finmarkiensis | Synoicum pulmonaria | Polycarpa fibrosa | Halocynthia pyriformis | Boltenia echinata', # Combine, hard to distinguish in the field. Ascidiacea is coarse, but represent a widely distributed group that would be vastly underestimated if left out. We therefore combine them at the 'Ascidiacea' level, only leaving 'Boltenia ovifera' and 'Eudistoma vitreum' as species, as they are easily identified
                     'Ophiura sp. | Ophiura sarsii | Ophiura robusta', # species easily distinguishable, but often found in similar locations and O. robusta is smaller and often hard to distinguish in trawl among important biomass of O. sarsii
                     'Lithodes maja | Lithodes sp.', # Very likely the same species
                     'Actinauge sp. | Actinauge cristata', # Very likely the same species, in any case will represent the same group. Combine at the genus level, as other less frequent species might be included. Keep in mind that the most likely species is the one being aggregated
                     'Actinostola sp. | Actinostola callosa', # Very likely the same species, in any case will represent the same group. Combine at the genus level, as other less frequent species might be included. Keep in mind that the most likely species is the one being aggregated
                     'Bolocera sp. | Bolocera tuediae', # Very likely the same species, in any case will represent the same group. Combine at the genus level, as other less frequent species might be included. Keep in mind that the most likely species is the one being aggregated
                     'Stephanauge sp. | Stephanauge nexilis', # Very likely the same species, in any case will represent the same group. Combine at the genus level, as other less frequent species might be included. Keep in mind that the most likely species is the one being aggregated
                     'Alcyonidium sp. | Alcyonidium pachydermatum', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Nuculana sp. | Nuculana tenuisulcata', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Mytilus sp. | Mytilus edulis', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Astarte sp. | Astarte subaequilatera | Astarte borealis', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Naticidae | Cryptonatica affinis | Euspira sp. | Euspira pallida', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements.
                     'Buccinum sp. | Buccinum scalariforme | Buccinum undatum', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Colus sp. | Colus stimpsoni | Colus pubescens | Plicifusus kroeyeri', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Neptunea sp. | Neptunea decemcostata | Neptunea despecta', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Ampelisca sp. | Ampelisca eschrichtii', # Very likely the same species, in any case will represent the same group. Combine at the genus level, as other less frequent species might be included.
                     'Dendronotus sp. | Dendronotus frondosus', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Tonicella sp. | Tonicella rubra', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Bathypolypus bairdii | Bathypolypus sp.', # Very likely the same species
                     'Nymphon sp. | Nymphon hirtipes', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Balanidae | Balanus balanus', # Keep even though coarse. The family represents the same type of species functionally and are hardly distinguishable in the field. Combined with Balanus balanus, which is often grouped in the Balanidae family during identification.
                     'Eualus gaimardii | Eualus gaimardii belcheri | Eualus gaimardii gaimardii', # Subspecies too detailed for my analyses
                     'Pagurus sp. | Pagurus pubescens', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Molpadia | Molpadia oolitica', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Gorgonocephalus sp. | Gorgonocephalus arcticus', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Eumicrotremus spinosus | Eumicrotremus spinosus variabilis', # Subspecies too detailed for my analyses
                     'Flabellum sp. | Flabellum alabastrum', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Epizoanthus sp. | Epizoanthus erdmanni', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Heliometra glacialis | Crinoidea', # Grouped at the species level. I is unlikely that there are other species in the St. Lawrence with the same physiology. If there are, it could be 'Poliometra proliza', which is found in the Arctic, but not recorded in the St. Lawrence. In any event, they share similar functions and are often found together in the Arctic.
                     'Margarites sp. | Margarites groenlandicus | Margarites costalis', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Boreotrophon sp. | Boreotrophon clathratus | Scabrotrophon fabricii', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Velutinidae | Velutina velutina | Limneria undata', # Combine, hard to distinguish in the field and likely composed of functionally similar species limited to the two listed here.
                     'Arrhoges occidentalis | Aporrhais sp.', # Combine at the species level, likely the only species of that group in the St. Lawrence
                     'Poraniomorpha sp. | Poraniomorpha hispida', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Musculus sp. | Musculus niger', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Cuspidaria sp. | Cuspidaria glacialis', # Identification: this group contains multiple species that generally have similar functional roles and habitat requirements, with C. glacialis likely included in Cuspidaria sp. Grouped at the genus level
                     'Gonostomatidae | Cyclothone microdon', # Combine at the family level. The only other species found in the St. Lawrence from this group is 'Cyclothone braueri', which has similar feeding strategies, habitat requirements and functional role.
                     'Myctophidae | Benthosema glaciale | Lampadena speculigera | Notoscopelus sp. | Notoscopelus elongatus', # Diverse group that can be hard to distinguish in the field and that are found in similar types of habitats. Most specimen are identified at the family level. Identifications at the genus or species level should be correct, but since there are vastly more speciments attributed to the family, all specimens from this family will be aggregated to the family level 'Myctophidae' and used for analyses. Maybe consider adding Neoscopelus macrolepidotus, which is also species found at great depths, but member of the family 'Neoscopelidae'. For now, consider individually.
                     'Gaidropsarus sp. | Gaidropsarus argentatus | Gaidropsarus ensis', # Identification: can be hard to distinguish species from this group and they are species that generally have similar functional roles and habitat requirements. G. ensis is generally located in the St. Lawrence, while G. argentatus is generally in the north-eastern Atlantic. Grouped at the genus level
                     'Artediellus sp. | Artediellus atlanticus | Artediellus uncinatus', # Identification: can be hard to distinguish species from this group in the field, except for mature males, and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Icelus sp. | Icelus bicornis | Icelus spatula', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Triglops sp. | Triglops murrayi | Triglops nybelini', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Myoxocephalus sp. | Myoxocephalus octodecemspinosus | Myoxocephalus scorpius', # Identification: can be hard to distinguish species from this group in the field and they are species that generally have similar functional roles and habitat requirements. Grouped at the genus level
                     'Lycodes terraenovae | Lycodes esmarkii', # Often confounded and might in fact be a single species. Grouped under the 'Lycodes terraenovae' umbrella, but keep in mind that it includes both species
                     'Ammodytes sp. | Ammodytes americanus | Ammodytes dubius' # Distribution and physiological overlap between the two species, record at the genus level.
                    )

    combineTaxa <- str_split(combineTaxa, ' \\| ')

    # Duplicate taxa names in northPluri to maatch taxon name chosen for aggregation
        for(i in 1:length(combineTaxa)) {
            nTaxa <- length(combineTaxa[[i]])
            # Identify records for taxon used for combination
                taxonCombRec <- northPluri[, 'N_EspSci'] == combineTaxa[[i]][1]

            # if taxon exist in list, replace taxa name that will be aggregated with taxonComb
            if(sum(taxonCombRec) > 0) {
                # Extract variables for taxon used for combination
                    taxonComb <- northPluri[taxonCombRec, ][1, c('EspGen','N_EspSci','N_EspF')]

                # Identify which taxa have to be combined with taxonComb
                    taxaReplace <- northPluri[, 'N_EspSci'] %in% combineTaxa[[i]][2:nTaxa]

                # Replace names of taxa to combine with taxonComb
                    northPluri[taxaReplace, c('EspGen','N_EspSci','N_EspF')] <- taxonComb
            } else {
                # Create new EspGen for taxon and create vector of variables used to describe taxon used for combination
                # !!!!! There may already be an existing number. Need to figure this out at some point !!!!!!
                    newEspGen <- max(as.numeric(paste(northPluri[,'EspGen']))) + 1
                    taxonComb <- data.frame(EspGen = newEspGen, N_EspSci = combineTaxa[[i]][1], N_EspF = combineTaxa[[i]][1], stringsAsFactors = FALSE)

                # Identify which taxa have to be combined with taxonComb
                    taxaReplace <- northPluri[, 'N_EspSci'] %in% combineTaxa[[i]][2:nTaxa]

                # Replace names of taxa to combine with taxonComb
                    northPluri[taxaReplace, c('EspGen','N_EspSci','N_EspF')] <- taxonComb
            }
        }

    # Aggregate data
        northPluri <- aggregate(cbind(WCapOri, SNb_Capt) ~ No_Rel + No_Stn + EspGen + N_EspSci + N_EspF + Source + Nbpc + DatDeTow + DatFiTow + EngGen + Resul + HreDeb + HreFin + LaDeTow + LoDeTow + LaFiTow + LoFiTow + Prof_1 + Prof_2, data = northPluri, FUN = sum, na.action = na.pass)

    # Order dataset by year, station and species number
        northPluri <- northPluri[order(northPluri[, 'No_Rel'], northPluri[, 'No_Stn'], northPluri[, 'EspGen']), ]

    # Individual ID
        ID <- seq(1,nrow(northPluri))
        northPluri <- cbind(ID, northPluri)

# List of taxa to remove
    removeTaxa <- c('Rhodophyta', # survey not designed to capture this taxon
                    'Laminaria sp.', # survey not designed to capture this taxon
                    'Cnidaria', # group is too coarse
                    'Anthozoa', # group is too coarse
                    'Nemertea', # group is too coarse & survey not designed to capture this taxon
                    'Priapulus caudatus', # survey not designed to capture this taxon
                    'Nematoda', # group is too coarse & survey not designed to capture this taxon
                    'Bivalvia', # group is too coarse
                    'Placopecten magellanicus', # coastal & shallow water species
                    'Scaphopoda', # group is too coarse
                    'Gastropoda', # group is too coarse
                    'Cephalopoda', # group is too coarse
                    'Onychoteuthidae', # group is too coarse, species are identified at the species level
                    'Lepidoteuthidae', # group is too coarse
                    'Polychaeta', # group is too coarse
                    'Maldanidae', # group is too coarse & survey not designed to capture this taxon
                    'Flabelligeridae', # group is too coarse
                    'Crustacea', # group is too coarse
                    'Cumacea', # group is too coarse
                    'Amphipoda', # group is too coarse
                    'Gammaridea', # group is too coarse
                    'Euphausiacea', # group is too coarse and survey not designed to capture this taxon
                    'Meganyctiphanes norvegica', # survey not designed to capture this taxon, pelagic species
                    'Mysida', # group is too coarse and survey not designed to capture this taxon
                    'Asteroidea', # group is too coarse
                    'Ophiuroidea', # group is too coarse
                    'Myctophiformes', # group is too coarse
                    'Gadidae', # group is too coarse
                    'Gadus sp.', # species differ in spatial distribution
                    'Hydrozoa', # group is too coarse and probably composed of pelagic hydrozoa P. lactea and S. martensii which can be confounded with scyphozoa A. aurita
                    'Scyphozoa', # group is too coarse and probably composed of pelagic scyphozoa A. aurita which can be confounded with hydrozoa P. lactea and S. martensii
                    'Brada inhabilis', # identification difficult, likely grouped as 'Polychaeta'
                    'Melinna cristata', # identification difficult, likely grouped as 'Polychaeta'
                    'Nereis pelagica', # identification difficult, likely grouped as 'Polychaeta'
                    'Nephtys sp.', # identification difficult, likely grouped as 'Polychaeta'
                    'Polyphysia crassa', # identification difficult, likely grouped as 'Polychaeta'
                    'Actiniaria', # group is too coarse
                    'Alcyonacea', # group is too coarse
                    'Ctenophora', # group is too coarse
                    'Turbellaria', # group is too coarse, survey not adapted
                    'Fecampiidae', # group is too coarse, survey not adapted
                    'Buccinidae', # goup is too coarse
                    'Nudibranchia', # group is too coarse
                    'Polyplacophora', # group is too coarse
                    'Sipuncula', # group is too coarse
                    'Echiura', # group it too coarse
                    'Pycnogonida', # group is too coarse
                    'Cirripedia', # group is too coarse
                    'Themisto compressa', # pelagic species, survey not adapted to properly capture
                    'Themisto libellula', # pelagic species, survey not adapted to properly capture
                    'Acanthephyra sp.', # pelagic species, survey not adapted to properly capture
                    'Acanthephyra pelagica', # pelagic species, survey not adapted to properly capture
                    'Eualus sp.', # species readily identified, hence too coarse at the genus level
                    'Lebbeus sp.', # species readily identified, hence too coarse at the genus level
                    'Spirontocaris sp.', # species readily identified, hence too coarse at the genus level
                    'Sabinea sp.', # species readily identified, hence too coarse at the genus level
                    'Gasterosteidae', # group is too coarse, alternatively, could be combined with Gasterosteus aculeatus aculeatus
                    'Holothuroidea', # group is too coarse
                    'Cephalaspidea', # group is too coarse
                    'Phaeophyceae', # survey not designed to capture this taxon
                    'Agarum sp.', # survey not designed to capture this taxon
                    'Agarum clathratum', # survey not designed to capture this taxon
                    'Nephtheidae', # group is too coarse, but ultimately could be grouped in a single group if species have similar distributions (Duva florida, Gersemia rubiformis, Drifa glomerata)
                    'Chondrus crispus', # survey not designed to capture this taxon
                    'Pteraster sp.', # group is too coarse, species easily distinguishable
                    'Asteriidae', # group is too coarse
                    'Scyliorhinidae', # group is too coarse
                    'Sternoptychidae', # group is too coarse
                    'Cottidae', # group is too coarse
                    'Liparidae', # group is too coarse
                    'Liparis sp.', # groups is too coarse, species distinguishable
                    'Paraliparis sp.', # groups is too coarse, species distinguishable
                    'Zoarcidae', # group is too coarse
                    'Lycenchelys sp.', # group is too coarse, species distinguishable
                    'Lycodes sp.', # group is too coarse, species distinguishable
                    'Stichaeidae', # group is too coarse
                    'Dentaliidae', # group is too coarse
                    'Isididae', # group is too coarse
                    'Cadlina laevis' # coastal species
                    )

    # spNames[!spNames[, 'N_EspSci'] %in% removeTaxa, ] # To visualize data removed

    # Remove taxa from northPluri
        northPluri <- northPluri[!northPluri[, 'N_EspSci'] %in% removeTaxa, ]

# Select taxa with sufficient data per year to be able to create models
    # What is the amount of data necessary?
    # Look at SDM literature, JSDM articles too

    # Count number of observations per species
        countRec <- numeric(nrow(northPluri))
        spNames <- aggregate(countRec ~ EspGen + N_EspSci + N_EspF, data = cbind(northPluri, countRec), FUN = length) %>%
                        .[order(.[, 'EspGen']), ]

    # Minimal number of records is 50 per year?
        minRec <- 50

    # Select species with nRec >= minRec
        spNames <- spNames[spNames[, 'countRec'] >= minRec, ]

    # Create subset of original northPluri dataset corrected for relevant taxa and number of records
        northPluriCor <- northPluri[northPluri[, 'N_EspSci'] %in% spNames[, 'N_EspSci'], ]

# Save analysis dataset
    saveRDS(northPluriCor, file = './RData/northPluriCor.rds')

# Visualize data informations
    # Structure of formatted dataset
        str(northPluriCor)

    # List of taxa with number of records
        spNames

    # Number of taxa
        nrow(spNames)

    # Average number of records and standard deviation
        mean(spNames[, 'countRec'])
        sd(spNames[, 'countRec'])

    # List number of stations
        numStations <- aggregate(No_Stn ~ No_Rel, data = unique(northPluriCor[, c('No_Rel','No_Stn')]), FUN = length)
        numStations

    # Total number of stations & Average number and standard deviation of stations per year
        sum(numStations[, 'No_Stn'])

        mean(numStations[, 'No_Stn'])

        sd(numStations[, 'No_Stn'])
