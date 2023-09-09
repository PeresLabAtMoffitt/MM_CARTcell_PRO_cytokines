# Import packages
library(tidyverse)


############################################################# Load data
path <- fs::path("", "Volumes", "Peres_Research","Myeloma", "PROs and cytokines")
path1 <- fs::path("", "Volumes", "Peres_Research","Myeloma", "Grant applications", 
                 "InNovation proposal_2022")


clinical <- 
  read_rds("~/Documents/GitHub/Peres/MM_CARTcell_mIF_exhaustion/clinical_with_new_variables.rds")

# discharge <- 
#   read_csv("/Volumes/Peres_Research/Myeloma/Analyses/ASH cytokines/Date of discharge_06062022.csv") %>% 
#   janitor::clean_names()
# discharge <- 
#   read_csv(paste0(path1, "/data/raw data/phi_cytokine_pt_data_07-06-2023.csv")) %>% 
#   janitor::clean_names() %>% 
#   select(mrn, initial_discharge_date)

cytokines <- 
  read.csv(paste0(path, "/data/Sept2023 R01 Prelim data",
                           "/ABECMA ELLA-Raw Data_LCP.csv")) %>% 
  janitor::clean_names()

pros <- 
  read.csv(paste0(path, "/data/Sept2023 R01 Prelim data",
                  "/pros_08082023.csv")) %>% 
  janitor::clean_names()


############################################################# Clean data
clinical <- clinical %>% 
  # left_join(., discharge %>% 
  #             mutate(mrn = as.character(mrn)),
  #           by = "mrn") %>% 
  # left_join(., discharge2 %>% 
  #             mutate(mrn = as.character(mrn)),
  #           by = "mrn") %>% 
  mutate_at(c("initial_discharge_date",
              "date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion"), 
            ~ as.Date(., format = "%m/%d/%y")) %>% 
  mutate(hospitalstay = initial_discharge_date - 
           date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion) %>% 
  # mutate(baseline_ferritin = case_when(
  #   baseline_ferritin >= 400                                    ~ "≥ ULN at LD",
  #   baseline_ferritin < 400                                    ~ "Normal"
  # )) %>% 
  mutate(peak_ferritin = case_when(
    max_ferritin >= 400                                    ~ "≥ ULN at LD",
    max_ferritin < 400                                    ~ "Normal"
  )) %>% 
  # mutate(baseline_crp = case_when(
  #   baseline_crp >= 0.5                                    ~ "≥ ULN at LD",
  #   baseline_crp < 0.5                                    ~ "Normal"
  # )) %>% 
  mutate(peak_crp = case_when(
    max_crp >= 0.5                                    ~ "≥ ULN at LD",
    max_crp < 0.5                                    ~ "Normal"
  ))

pros <- pros %>% 
  mutate(mrn = as.character(mrn)) %>% 
  mutate(day = case_when(
    day == -1          ~ "baseline",
    TRUE               ~ as.character(day)
  )) %>% 
  pivot_wider(names_from = day,
              values_from = -c(mrn, day),
              names_glue = "{.value}_{day}"
  )
pros <- pros %>%
  mutate(across(contains("factg_total"), ~ case_when(
    . < 63                                              ~ "Low",
    . >= 63                                             ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("factg_pwb"), ~ case_when(
    . < 16                                              ~ "Low",
    . >= 16                                             ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("factg_swb"), ~ case_when(
    . < 17                                              ~ "Low",
    . >= 17                                             ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("factg_ewb"), ~ case_when(
    . < 14                                              ~ "Low",
    . >= 14                                             ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("factg_fwb"), ~ case_when(
    . < 12                                              ~ "Low",
    . >= 12                                             ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("ti_"), ~ case_when(
    . == 0                                              ~ "None",
    . < 2                                               ~ "Mild",
    . < 3                                               ~ "Moderate",
    . < 4                                               ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>% 
  mutate(across(contains("ti_"), ~ case_when(
    . < 2                                               ~ "None-Mild",
    . < 3                                               ~ "Moderate",
    . < 4                                               ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat2")) %>% 
  mutate(across(contains("ti_"), ~ case_when(
    . < 3                                               ~ "< 3",
    . >= 3                                              ~ "≥ 3",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat3")) %>% 
  mutate(across(contains("depression_tscore"), ~ case_when(
    . < 55                                              ~ "Normal",
    . <= 59                                             ~ "Mild",
    . <= 69                                             ~ "Moderate",
    . >= 70                                             ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(promis_depression_tscore_baseline_cat = 
           factor(promis_depression_tscore_baseline_cat, 
                  levels = c("Normal", "Mild",
                             "Moderate"))
  ) %>% 
  mutate(across(contains("anxiety_tscore"), ~ case_when(
    . < 55                                              ~ "Normal",
    . <= 59                                             ~ "Mild",
    . <= 69                                             ~ "Moderate",
    . >= 70                                             ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(promis_anxiety_tscore_baseline_cat = 
           factor(promis_anxiety_tscore_baseline_cat, 
                  levels = c("Normal", "Mild",
                             "Moderate"))
  ) %>% 
  mutate(across(contains("fatigue_tscore"), ~ case_when(
    . < 55                                              ~ "Normal",
    . <= 59                                   ~ "Mild",
    . <= 69   ~ "Moderate",
    . >= 70   ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("fatigue_tscore"), ~ case_when(
    . < 55                                              ~ "Normal",
    . <= 59                                   ~ "Mild",
    . > 59   ~ "Moderate-Severe",
    # . >= 70   ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat2")) %>%
  mutate(across(contains("fatigue_tscore"), ~ case_when(
    . < 60                                    ~ "< 60",
    . >= 60                                   ~ "≥ 60",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat3")) %>%
  mutate(across(contains("pain_tscore"), ~ case_when(
    . < 55                                              ~ "Normal",
    . <= 59                                   ~ "Mild",
    . <= 69   ~ "Moderate",
    . >= 70   ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("pain_tscore"), ~ case_when(
    . < 55                                              ~ "Normal",
    . <= 59                                   ~ "Mild",
    . > 59   ~ "Moderate-Severe",
    # . >= 70   ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat2")) %>%
  mutate(across(contains("pain_tscore"), ~ case_when(
    . < 39                                              ~ "< 39",
    . >= 39                                             ~ "≥ 39",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat3")) %>%
  mutate(across(contains("sleep_tscore"), ~ case_when(
    . < 55                                              ~ "Normal",
    . <= 59                                   ~ "Mild",
    . <= 69   ~ "Moderate",
    . >= 70   ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("physical_tscore"), ~ case_when(
    . < 30                                              ~ "Severe",
    . <= 39                                   ~ "Moderate",
    . <= 44   ~ "Mild",
    . >= 45   ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("physical_tscore"), ~ case_when(
    # . < 30                                              ~ "Severe",
    . <= 39                                             ~ "Moderate-Severe",
    . <= 44                                             ~ "Mild",
    . >= 45                                             ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat2")) %>%
  mutate(across(contains("physical_tscore"), ~ case_when(
    . < 39                                              ~ "< 39",
    . >= 39                                             ~ "≥ 39",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat3")) %>%
  mutate(across(contains("social_tscore"), ~ case_when(
    . < 30                                              ~ "Severe",
    . <= 39                                   ~ "Moderate",
    . <= 44   ~ "Mild",
    . >= 45   ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(promis_social_tscore_baseline_cat = 
           factor(promis_social_tscore_baseline_cat, 
                  levels = c("Normal", "Mild",
                             "Moderate", "Severe"))
  ) %>% 
  mutate(across(contains("cognitive_tscore"), ~ case_when(
    . < 30                                              ~ "Severe",
    . <= 39                                   ~ "Moderate",
    . <= 44   ~ "Mild",
    . >= 45   ~ "Normal",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) %>%
  mutate(across(contains("globalpain"), ~ case_when(
    . == 0                                              ~ "None",
    . <= 4                                   ~ "Mild",
    . <= 7   ~ "Moderate",
    . >= 8   ~ "Severe",
    TRUE                                                ~ NA_character_)
    , .names = "{.col}_cat")) 


cytokines <- cytokines %>% 
  mutate(mrn = as.character(mrn)) %>% 
# Determining baseline/day-6 cytokine values
  mutate(IL2_baseline = case_when(
    timepoint == "Day-6"            ~ il2
  )) %>% 
  mutate(IL6_baseline = case_when(
    timepoint == "Day-6"            ~ il6
  ), .after = il6) %>% 
  mutate(IFNg_baseline = case_when(
    timepoint == "Day-6"            ~ if_ng
  )) %>% 
  mutate(TNFa_baseline = case_when(
    timepoint == "Day-6"            ~ tn_fa
  )) %>% 
  inner_join(clinical %>% 
              select(mrn, hospitalstay), 
             .,
            by = "mrn") %>% 
  group_by(mrn) %>% 
  mutate(peak_il6 = case_when(
    days_of_car_t > 0 & 
      days_of_car_t <= hospitalstay      ~ max(il6, na.rm = TRUE)
  ), .after = il6) %>% 
  mutate(peak_il2 = case_when(
    days_of_car_t > 0 & 
      days_of_car_t <= hospitalstay      ~ max(il2, na.rm = TRUE)
  ), .after = il2) %>% 
  mutate(peak_ifng = case_when(
    days_of_car_t > 0 & 
      days_of_car_t <= hospitalstay      ~ max(if_ng, na.rm = TRUE)
  ), .after = if_ng) %>% 
  mutate(peak_tnfa = case_when(
    days_of_car_t > 0 & 
      days_of_car_t <= hospitalstay      ~ max(tn_fa, na.rm = TRUE)
  ), .after = tn_fa) %>% 
  fill(c(starts_with("peak"), ends_with("baseline")), .direction = "updown") %>% 
  ungroup() %>% 
  select(mrn, starts_with("peak"), ends_with("baseline")) %>% 
  distinct()


############################################################# Merge data
data <- full_join(clinical %>% 
                    mutate(data_has_clinical = "Yes"), 
                  pros %>% 
                    mutate(data_has_pros = "Yes"),
                  by = "mrn") %>% 
  full_join(cytokines %>% 
              mutate(data_has_cytokines = "Yes"),
            by = "mrn") %>% 
  select(mrn, starts_with("data_has_"), everything())

write_rds(data, "pros_cytokines_data.rds")


# End cleaning

