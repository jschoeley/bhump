# Derive microdata in event-history format for analysis of
# fetal-infant transition

# Init --------------------------------------------------------------------

library(qs2)
library(readr)
library(data.table)
library(lubridate)
library(magrittr)

# use all available CPUs
setDTthreads(0)

paths <- list()
paths$input <- list(
  fetoinfant.qs = 'tmp/20-fetoinfant.qs',
  cod.csv = 'dat/10-cod-list/cod.csv'
)
paths$output <- list(
  fetoinfant.qs = 'tmp/21-fetoinfant.qs'
)

# constants
cnst <-
  list(
    left_truncation_gestage = 24,
    right_censoring_gestage = 76.99
  )

# Load data ---------------------------------------------------------------

# US fetal- infant deaths and births individual level data
fetoinfant <- qs_read(paths$input$fetoinfant.qs)
setDT(fetoinfant)

# Subset ------------------------------------------------------------------

# only consider cases which survived to
# week 24, i.e. were delivered (alive or dead)
# at week 24 or later and were conceived during
# the years 1989, 1999, 2009, and 2014.
# to that end we only consider deliveries during
# 1989 & 1990 -> 1989,
# 1999 & 2000 -> 1999,
# 2009 & 2010 -> 2009
# 2014 & 2015 -> 2014
fetoinfant <-
  fetoinfant[
    gestation_at_delivery_w >= cnst$left_truncation_gestage &
      date_of_delivery_y %in%
      c(1989, 1990,
        1999, 2000,
        2009, 2010,
        2014, 2015),
    ]

#fetoinfant[, .(count = .N), by = date_of_delivery_y]

# Calculate date of conception --------------------------------------------

# calculate date of conception as
# date of delivery - (weeks of gestation at delivery - 2 weeks)
# the second term is the ferilization age at delivery which
# is shorter than the gestational age as fertilization on average
# happens 2 weeks following the last period
fetoinfant[
  ,
  date_of_delivery_ym :=
    paste(date_of_delivery_y,
          date_of_delivery_m,
          '15', sep = '-') |>
    ymd()
  ]
fetoinfant[
  ,
  date_of_conception_ym :=
    (date_of_delivery_ym - weeks(gestation_at_delivery_w - 2)) |>
    round_date(unit = 'month')
  ]
fetoinfant[
  ,
  date_of_conception_y :=
    year(date_of_conception_ym)
  ]

# subset to conception cohorts 1989, 1999, 2009, and 2014
fetoinfant <-
  fetoinfant[
    date_of_conception_y %in% c(1989, 1999, 2009, 2014)
  ]

# delete space after icd code
fetoinfant[, cod_icd10 := gsub(' ', '', cod_icd10)]

#fetoinfant[, .(count = .N), by = date_of_conception_y]
dcast(
  fetoinfant[, .(count = .N),
             by = .(date_of_conception_y, date_of_delivery_y, type)],
  date_of_conception_y + date_of_delivery_y ~ type,
  value.var = 'count'
)

# Recode cause of death categories ----------------------------------------

# load mapping between cod categories and icd-10 codes
cod_codes <- read_csv(paths$input$cod.csv)
cod_codes <- as.data.table(cod_codes)
# select relevant cod coding here
cod_codes <- cod_codes[,.(cod_icd10, cod_cat = cod_cat)]

# add cod categories to data
fetoinfant <-
  merge(fetoinfant, cod_codes, by = 'cod_icd10', all.x = TRUE)

# define special mappings
fetoinfant[
  ,
  cod_cat := fcase(
    # prematurity only if gestation at birth was early
    (startsWith(cod_icd10, c('P22', 'P26', 'P27', 'P28')) &
       gestation_at_delivery_w >= 37), 'otherspecific',
    # code explicit COD NAs for fetal death
    (type == 'fetus' & (cod_icd10 == '' | is.na(cod_icd10))), 'unknown',
    # code explicit COD NAs for infant deaths
    (type == 'infant' & !is.na(age_at_death_d) & (cod_icd10 == '' | is.na(cod_icd10))), 'unknown',
    # defaults
    !(cod_icd10 == '' | is.na(cod_icd10)), cod_cat,
    default = NA
  )
]

# Add flags for vital events ----------------------------------------------

# add flags for vital events fetal-death, life-birth,
# neonatal and post-neonatal death and survival
fetoinfant[
  ,
  `:=`(
    # fetal deaths are all cases from the
    # fetal death file
    fetal_death =
      type == 'fetus',
    # life births are all cases from the
    # birth registry file
    life_birth =
      type == 'infant',
    # neonatal deaths are all cases with an entry
    # for the age at death in days which is < 7
    # FALSE & NA = FALSE
    neonatal_death =
      (!is.na(age_at_death_d)) & (age_at_death_d < 7),
    # neonatal survivors are all cases from the infant file
    # without an age at death entry or with age at death >= 7
    neonatal_survivor =
      type == 'infant' & (is.na(age_at_death_d) | (age_at_death_d >= 7)),
    # postneonatal deaths are all cases with an
    # age at death >= 7
    postneonatal_death =
      (!is.na(age_at_death_d)) & (age_at_death_d >= 7),
    # postneonatal survivors are all cases from the infant file
    # without an age at death
    postneonatal_survivor =
      type == 'infant' & is.na(age_at_death_d)
  )
]

# Add timing of vital events ----------------------------------------------

# assume uniform distribution of events over week of
# gestation for fetal deaths and life-births
fetoinfant[
  ,
  `:=`(
    # gestational age at fetal death in weeks
    # assume uniform distribution of fetal deaths
    # within week
    gestage_at_fetal_death_w =
      fifelse(
        fetal_death,
        gestation_at_delivery_w + runif(.N),
        as.numeric(NA)
      ),
    # gestational age at life-birth
    # assume uniform distribution of life-births
    # within week
    gestage_at_life_birth_w =
      fifelse(
        life_birth,
        gestation_at_delivery_w + runif(.N),
        as.numeric(NA)
      ),
    # chronological age at infant death in (fractional) weeks
    # assume that death occurs uniformly throughout a day
    age_at_death_w =
      (age_at_death_d + runif(.N))/7
  )]

fetoinfant[
  ,
  `:=`(
    # gestational age at neonatal death in weeks
    gestage_at_neonatal_death_w =
      fifelse(
        neonatal_death,
        gestage_at_life_birth_w + age_at_death_w,
        as.numeric(NA)
      ),
    # gestational age at postneonatality in weeks
    gestage_at_postneonatality_w =
      fifelse(
        neonatal_survivor,
        gestage_at_life_birth_w + 1,
        as.numeric(NA)
      ),
    # gestation at postneonatal death in weeks
    gestage_at_postneonatal_death_w =
      fifelse(
        postneonatal_death,
        gestage_at_life_birth_w + age_at_death_w,
        as.numeric(NA)
      ),
    # censoring at end of week 76 = 1 year post viability
    gestage_at_censoring_w =
      cnst$right_censoring_gestage
  )]

# add id variable
fetoinfant[,id := 1:.N]

# select variables
fetoinfant <-
  fetoinfant[,
          .(
            id, sex,
            race_and_hispanic_orig_of_mother,
            education_of_mother,
            date_of_conception_y, date_of_conception_ym,
            fetal_death, life_birth,
            neonatal_death, neonatal_survivor,
            postneonatal_death, postneonatal_survivor,
            gestage_at_fetal_death_w,
            gestage_at_life_birth_w,
            gestage_at_neonatal_death_w,
            gestage_at_postneonatality_w,
            gestage_at_postneonatal_death_w,
            gestage_at_censoring_w,
            cod_cat
          )]

# Convert data to event-history format ------------------------------------

# population strata
strata <-
  c(
    'id', 'sex',
    'race_and_hispanic_orig_of_mother',
    'education_of_mother',
    'date_of_conception_y'
  )

# Fetal death event histories ---------------------------------------------

# fetal death event histories
# fetus -> death
fetal_deaths <- fetoinfant[fetal_death == TRUE]
fetal_deaths_histories <-
  cbind(
    fetal_deaths[,..strata],
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'death',
      destination2 = fetal_deaths$cod_cat,
      exit_time = fetal_deaths$gestage_at_fetal_death_w
    )
  )

# Neonatal death event histories ------------------------------------------

# neonatal death event histories
neonatal_deaths <- fetoinfant[neonatal_death == TRUE]
neonatal_deaths_histories <-
  rbind(
    # fetus -> neonate
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'neonate',
      destination2 = 'neonate',
      exit_time = neonatal_deaths$gestage_at_life_birth_w
    ),
    # neonate -> death
    data.table(
      origin = 'neonate',
      entry_time = neonatal_deaths$gestage_at_life_birth_w,
      destination = 'death',
      destination2 = neonatal_deaths$cod_cat,
      exit_time = neonatal_deaths$gestage_at_neonatal_death_w
    )
  )
# why does this work?
# the number of rows in neonatal_deaths_histories is an integer
# multiple of the number of rows in neonatal_deaths, the latter
# rows being recycled to the length of the neonatal_deaths_histories.
# as all the individuals in neonatal_deaths_histories experience the
# same number of events/rows and the order of events is encoded in the
# order of rows, the recycled id's actually belong to the same individual.
# the same logic applies to the other event_histories.
neonatal_deaths_histories <-
  cbind(
    neonatal_deaths[,..strata],
    neonatal_deaths_histories
  )

# Postneonatal death event histories --------------------------------------

# postneonatal death event histories
# fetus -> neonatal -> postneonatal -> (death|censored)
postneonatal_deaths <- fetoinfant[postneonatal_death == TRUE]
postneonatal_deaths[
  ,
  censored := gestage_at_postneonatal_death_w > gestage_at_censoring_w
  ]
postneonatal_deaths_histories <-
  rbind(
    # fetus -> neonate
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'neonate',
      destination2 = 'neonate',
      exit_time = postneonatal_deaths$gestage_at_life_birth_w
    ),
    # neonate -> post-neonate
    data.table(
      origin = 'neonate',
      entry_time = postneonatal_deaths$gestage_at_life_birth_w,
      destination = 'postneonate',
      destination2 = 'postneonate',
      exit_time = postneonatal_deaths$gestage_at_postneonatality_w
    ),
    # post-neonate -> (death | censoring)
    # infant -> (death|censoring)
    data.table(
      origin = 'postneonate',
      entry_time = postneonatal_deaths$gestage_at_postneonatality_w,
      destination =
        fifelse(
          postneonatal_deaths$censored == FALSE,
          'death',
          'censored'
        ),
      destination2 =
        fifelse(
          postneonatal_deaths$censored == FALSE,
          postneonatal_deaths$cod_cat,
          'censored'
        ),
      exit_time =
        fifelse(
          postneonatal_deaths$censored == FALSE,
          postneonatal_deaths$gestage_at_postneonatal_death_w,
          postneonatal_deaths$gestage_at_censoring_w
        )
    )
  )
postneonatal_deaths_histories <-
  cbind(
    postneonatal_deaths[,..strata],
    postneonatal_deaths_histories
  )

# Postneonatal survivor event histories -----------------------------------

# infant survivor event histories
# fetus -> neonatal -> postneonatal -> censored
postneonatal_survivors <- fetoinfant[postneonatal_survivor == TRUE]
postneonatal_survivors_histories <-
  rbind(
    # fetus -> neonate
    data.table(
      origin = 'fetus',
      entry_time = cnst$left_truncation_gestage,
      destination = 'neonate',
      destination2 = 'neonate',
      exit_time = postneonatal_survivors$gestage_at_life_birth_w
    ),
    # neonate -> post-neonate
    data.table(
      origin = 'neonate',
      entry_time = postneonatal_survivors$gestage_at_life_birth_w,
      destination = 'postneonate',
      destination2 = 'postneonate',
      exit_time = postneonatal_survivors$gestage_at_postneonatality_w
    ),
    # post-neonate -> (death | censoring)
    # infant -> (death|censoring)
    data.table(
      origin = 'postneonate',
      entry_time = postneonatal_survivors$gestage_at_postneonatality_w,
      destination = 'censored',
      destination2 = 'censored',
      exit_time = postneonatal_survivors$gestage_at_censoring_w
    )
  )
postneonatal_survivors_histories <-
  cbind(
    postneonatal_survivors[,..strata],
    postneonatal_survivors_histories
  )

# Complete feto-infant event histories ------------------------------------

# event histories for each subject during
# the feto-infant period
fetoinfant_event_histories <-
  rbind(
    fetal_deaths_histories,
    neonatal_deaths_histories,
    postneonatal_deaths_histories,
    postneonatal_survivors_histories
  )
fetoinfant_event_histories <-
  fetoinfant_event_histories[order(id)]

# Consistency checks ------------------------------------------------------

# cohort size matches
stopifnot(
  'Cohort sizes do not match between individual and event history formats' =
    fetoinfant[,.N] == fetoinfant_event_histories[,length(unique(id))],
  'Fetal deaths do not match between individual and event history formats' =
    fetoinfant_event_histories[
      origin == 'fetus' & destination == 'death',
      .N
    ] ==
    fetoinfant[
      fetal_death == TRUE,
      .N
    ],
  'Neonatal deaths do not match between individual and event history formats' =
    fetoinfant_event_histories[
      origin == 'neonate' & destination == 'death',
      .N
    ] ==
    fetoinfant[
      neonatal_death == TRUE,
      .N
    ],
  'Post-neonatal deaths do not match between individual and event history formats' =
    fetoinfant_event_histories[
      origin == 'postneonate' & destination == 'death',
      .N
    ] ==
    fetoinfant[
      postneonatal_death == TRUE &
        (gestage_at_postneonatal_death_w < gestage_at_censoring_w),
      .N
    ],
  'Censorings do not match between individual and event history formats' =
    fetoinfant_event_histories[
      origin == 'postneonate' & destination == 'censored',
      .N
    ] ==
    fetoinfant[
      postneonatal_survivor == TRUE |
        (gestage_at_postneonatal_death_w >= gestage_at_censoring_w),
      .N
    ],
  'Entry time is not always smaller than exit time' =
    all(fetoinfant_event_histories$entry_time <
          fetoinfant_event_histories$exit_time),
  'Spell exit and entry times do not line up' =
    fetoinfant_event_histories[
      .N > 1,
      .(same = entry_time[-1] == exit_time[-length(exit_time)]),
      by = id
    ][,
      all(same)
    ],
  'Observed CODs are not all featured in the cod table' =
    all(unique(fetoinfant$cod_icd10) %in% cod_codes$cod_icd10)
)

# Export ------------------------------------------------------------------

# save the processed microdata
qs_save(
  fetoinfant_event_histories,
  file = paths$output$fetoinfant.qs
)
