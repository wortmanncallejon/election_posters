###########################################

# PROJECT TITLE: Reinickendorf Wiederholungswahl Plakatanalyse

# AUTHOR: Felix Wortmann Callejón

# Date: 2024-02-12

###########################################

# Load Packages and Functions ----
setwd('/Users/callejon/Desktop/Privat/GitHub Repos/election_posters/')

pacman::p_load("here", "dplyr", "readr", "ggplot2", "ebal", "stringr", "fixest", "MatchIt", "sf", "purrr")

run_bern <- function(n_trials, n_draws, prob, ...) unlist(purrr::pmap(list(rep(n_draws, n_trials), rep(prob, n_trials)), function(d, p) mean(rbinom(d, 1, p))))
get_weights <- function(df, varlist) {
  
  covars = colnames(df)[varlist]
  X <- as.matrix(df[,covars])
  XX <- ebal::matrixmaker(X)
  
  
  # Construct vector of colinear covars
  cvs <- character(0)
  
  for (cv in covars) {
    cvs = c(cvs,paste0(cv,".",cv))
  }
  
  W <- ifelse(df$D == 1, T, F)
  
  X <- XX[, !colnames(XX) %in% cvs]
  
  out.eb <- ebal::ebalance(Treatment = W, X = X)
  
  tibble(id = df$id, 
         D = df$D) %>% 
    mutate(weight = ifelse(D == 1,1,out.eb$w)) %>% 
    select(-D) %>% 
    return()
}
signif <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", T ~ "")
calc_coefs <- function(df) {
  summarise(df,
            w_est = summary(lm(value ~ D, weights = weight))[["coefficients"]][2,1],
            w_se = summary(lm(value ~ D, weights = weight))[["coefficients"]][2,2],
            est = summary(lm(value ~ D))[["coefficients"]][2,1],
            se = summary(lm(value ~ D))[["coefficients"]][2,2],
            df = summary(lm(value ~ D))[["df"]][2],
            .by = name) %>% 
    return()
}

# Load Data -----

treated <- pull(distinct(read_csv(here("data", "data_rbb.csv")), Briefwahlbezirksnummer))

read_2021_data <- function() {
  
  zweit <- readxl::read_xlsx(here("data", "DL_BE_BU2021.xlsx"), sheet = 4) %>% 
    select(Adresse, `Bundestags-\r\nwahlkreis`, Briefwahlbezirk, `Wahlberechtigte insgesamt`, Wählende, `Gültige Stimmen`, CDU:Volt) %>%
    summarise(across(`Wahlberechtigte insgesamt`:Volt, \(x) sum(x, na.rm = T)), .by = c("Bundestags-\r\nwahlkreis", "Briefwahlbezirk")) %>% 
    mutate(wdh = ifelse(Briefwahlbezirk %in% treated,1,0)) %>% 
    tidyr::pivot_longer(cols = c(CDU:Volt),
                        names_to = "Partei", values_to = "Stimmen") %>% 
    janitor::clean_names() %>% 
    mutate(anteil = stimmen / gultige_stimmen,
           partei = ifelse(partei %in% c( "Die PARTEI", "Tierschutzpartei", "PIRATEN", "Die Grauen", "FREIE WÄHLER", "Gesundheitsforschung", "ÖDP", "du.", "V-Partei³", "DKP", "MLPD", "BüSo", "SGP", "LKR", "NPD", "Die Humanisten", "Team Todenhöfer", "Volt"), "Sonstige", partei)) %>% 
    summarise(across(wahlberechtigte_insgesamt:wdh, \(x) mean(x, na.rm = T)),
              across(stimmen:anteil, \(x) sum(x, na.rm =T)),
              .by = c("bundestags_wahlkreis", "briefwahlbezirk", "partei")) %>% 
    mutate(stimme = "Zweitstimme") %>% 
    filter(wdh == 1) 
  
  erst <- readxl::read_xlsx(here("data", "DL_BE_BU2021.xlsx"), sheet = 3) %>% 
    select(Adresse, `Bundestags-\r\nwahlkreis`, Briefwahlbezirk, `Wahlberechtigte insgesamt`, Wählende, `Gültige Stimmen`, CDU:`Anderer K....56`) %>%
    summarise(across(`Wahlberechtigte insgesamt`:`Anderer K....56`, \(x) sum(x, na.rm = T)), .by = c("Bundestags-\r\nwahlkreis", "Briefwahlbezirk")) %>% 
    mutate(wdh = ifelse(Briefwahlbezirk %in% treated,1,0)) %>% 
    tidyr::pivot_longer(cols = c(CDU:`Anderer K....56`),
                        names_to = "Partei", values_to = "Stimmen") %>% 
    janitor::clean_names() %>% 
    mutate(anteil = stimmen / gultige_stimmen,
           partei = ifelse(partei %in% c("SPD", "CDU", "GRÜNE", "FDP", "AfD", "DIE LINKE"), partei, "Sonstige")) %>% 
    summarise(across(wahlberechtigte_insgesamt:wdh, \(x) mean(x, na.rm = T)),
              across(stimmen:anteil, \(x) sum(x, na.rm =T)),
              .by = c("bundestags_wahlkreis", "briefwahlbezirk", "partei")) %>% 
    mutate(stimme = "Erststimme") %>% 
    filter(wdh == 1) 
  
  return(bind_rows(erst,zweit))
  
}
read_2024_data <- function() {
  zweit <- read_delim(here("data","Datenexport_BUNDESTAGSWAHL2024_Zweitstimme_W_BE.csv"), delim = ";") %>% 
    select(-matches("P\\d{2}p")) %>% 
    rename(CDU = P01,
           "DIE LINKE" = P02,
           SPD = P03,
           GRÜNE = P04,
           AfD = P05,
           FDP = P06,
           id = Briefwahlbezirk,
           gueltige_stimmen = Gueltig) %>%
    mutate(Sonstige = rowSums(select(., P07:P60))) %>% 
    select(id, gueltige_stimmen, CDU:FDP, Sonstige) %>% 
    summarise(across(gueltige_stimmen:Sonstige, \(x) sum(x, na.rm = T)), .by = id) %>% 
    tidyr::pivot_longer(cols = c(CDU:Sonstige), names_to = "partei", values_to = "sim") %>% 
    mutate(sim = sim/gueltige_stimmen,
           stimme = "Zweitstimme") %>% 
    filter(id %in% treated) 
  
  erst <- read_delim(here("data","Datenexport_BUNDESTAGSWAHL2024_Erststimme_W_BE.csv"), delim = ";") %>% 
    select(-matches("P\\d{2}p")) %>% 
    rename(CDU = P01,
           "DIE LINKE" = P02,
           SPD = P03,
           GRÜNE = P04,
           AfD = P05,
           FDP = P06,
           id = Briefwahlbezirk,
           gueltige_stimmen = Gueltig) %>% 
    mutate(Sonstige = rowSums(select(., P07:P60))) %>% 
    select(id, gueltige_stimmen, CDU:FDP, Sonstige) %>% 
    summarise(across(gueltige_stimmen:Sonstige, \(x) sum(x, na.rm = T)), .by = id) %>% 
    tidyr::pivot_longer(cols = c(CDU:Sonstige), names_to = "partei", values_to = "sim") %>% 
    mutate(sim = sim/gueltige_stimmen,
           stimme = "Erststimme") %>% 
    filter(id %in% treated) 
  
  return(bind_rows(erst,zweit))
}

data <- read_2021_data()
data_2024 <- read_2024_data()

social <- readxl::read_xls(here("data", "Sozialstrukturdaten_Wiederholung.xls"), sheet = 1) %>% 
  filter(Brief %in% treated) %>% 
  mutate(D = factor(ifelse(Bezirksname == "Reinickendorf",1,0))) %>% 
  rename(id = Brief) %>% 
  select(id, D, everything(),  -Bezirksname:-`_merge`)

bal_vars <- c(4:5, 7:8)

wdh_data <- data %>% 
  #filter(stimme == "Zweitstimme") %>% 
  mutate(D = factor(ifelse(bundestags_wahlkreis == "77",1,0), levels = c(0,1))) %>% 
  inner_join(get_weights(social, bal_vars), by = c("briefwahlbezirk" = "id")) %>% 
  select(briefwahlbezirk, D, weight, gultige_stimmen, partei, anteil, stimme)
  


wdh_data <- wdh_data %>% 
  left_join(select(data_2024, id, partei, sim, stimme), by = c("briefwahlbezirk" = "id", "partei", "stimme"))


# Matching ----

bal_data <- social %>% 
  inner_join(get_weights(social, bal_vars), by = "id") %>% 
  tidyr::pivot_longer(-c(id,D,weight)) %>% 
  mutate(value = value/100) %>% 
  calc_coefs() %>% 
  rename(Variables = name) %>%  
  tidyr::pivot_longer(-c(Variables,df)) %>% 
  mutate(Matched = ifelse(str_detect(name, "w_"),"Gewichtet", "Ungewichtet"),
         Statistic = ifelse(str_detect(name, "est"), "est", "se")) %>% 
  dplyr::select(-name) %>% 
  tidyr::pivot_wider(names_from = "Statistic", values_from = "value") %>% 
  mutate(t_crit = qt(0.025, df, lower.tail = F),
         lwr = est - t_crit*se,
         upr = est + t_crit*se,
         in_eb = ifelse(Variables %in% colnames(social)[bal_vars],"Ja","Nein"))

top_diffs <- bal_data %>% 
  filter(Matched == "Ungewichtet") %>% 
  arrange(desc(abs(est))) %>% 
  slice(1:12) %>% 
  select(Variables) %>% pull() %>% paste(collapse = " + ")

m.out <- matchit(as.formula(paste0("D ~ ", top_diffs)), data = social, method = "optimal")
data_matched <- match.data(m.out)

wdh_data %>% 
  filter(briefwahlbezirk %in% data_matched$id) %>% 
  summarise(est = summary(lm(sim ~ D))[["coefficients"]][2,1],
            se = summary(lm(sim ~ D))[["coefficients"]][2,2],
            df = summary(lm(sim ~ D))[["df"]][2],
            .by = c("partei", "stimme")) %>%
  mutate(t_crit = qt(0.025, df, lower.tail = F),
         lwr = est - t_crit*se,
         upr = est + t_crit*se,
         partei = factor(partei, levels = rev(c("SPD", "CDU", "GRÜNE", "AfD", "FDP", "DIE LINKE", "Sonstige")))) %>%
  ggplot(aes(est, partei,xmin = lwr, xmax = upr, color = partei)) +
  scale_x_continuous("ATT Schätzungen", labels = function(x) paste0(round(x*100,0),"%P")) +
  scale_y_discrete(NULL) +
  scale_color_manual(values = rev(c("#CC0000", "#000000", "#0E8C1D", "#005EA4", "#FFC000", "#CC0066", "#BFBFBF")), guide = "none") +
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_pointrange(size = 1, linewidth = 1) +
  facet_wrap(~stimme) +
  theme_light(base_size = 18) +
  theme(panel.grid = element_blank()) +
  ggtitle('Matching-Schätzung des "Plakateffekts"')


ggsave(file = "matching.png", width = 25, height = 15, units = "cm", dpi = 600)

# DiD -----
get_did_counterfactual <- function(s) {
  
  fixest.model <- feols(anteil ~ D*`T`, data = filter(did_data, partei == "GRÜNE" & stimme == s)) 
  coefs <- fixest.model[["coeftable"]]
  
  tibble(`T` = c(21,21,24,24,24,21),
         D = factor(c(0,1,0,1,1,1)),
         Counterfactual = c(0,0,0,0,1,1)) %>% 
    mutate(Y = ifelse(D == 0, coefs[1,1], coefs[1,1] + coefs[2,1]) + ifelse(`T` == 21, 0, coefs[3,1]) + ifelse(Counterfactual == 1 & `T` == 24, coefs[4,1], 0),
           D = factor(case_when(D == 0 ~ "Kontrollgruppe",
                                D == 1 & Counterfactual == 0 ~ "Reinickendorf (Kontrafaktisch)",
                                D == 1 & Counterfactual == 1 ~ "Reinickendorf (Beobachtet)"),
                      levels = c("Kontrollgruppe", "Reinickendorf (Beobachtet)", "Reinickendorf (Kontrafaktisch)")),
           stimme = s) %>% 
    return()
}
get_text_data <- function(s) {
  model <- feols(anteil ~ D*`T`, data = filter(did_data, partei == "GRÜNE" & stimme == s), cluster = "briefwahlbezirk") 
  
  tibble(x = c(24.1),
         y = c(model[["coeftable"]][4,1]/2 + sum(model[["coeftable"]][1:3,1])),
         lab = c(paste0(round(model[["coeftable"]][4,1]*100,2),"%P ",signif(model[["coeftable"]][4,4]))),
         stimme = s) %>% 
    return()
}
get_brace_data <- function(s) {
  model <- feols(anteil ~ D*`T`, data = filter(did_data, partei == "GRÜNE" & stimme == s)) 
  
  tibble(x = c(24.001,24.1),
         y = c(sum(model[["coeftable"]][1:3,1]), sum(model[["coeftable"]][1:4,1])),
         stimme = s) %>% 
    return()
}

did_data <- data_2024 %>% 
  select(id, partei, stimme, sim) %>% 
  rename(briefwahlbezirk = id,
         anteil = sim) %>% 
bind_rows(select(data, briefwahlbezirk, partei, stimme, anteil),
          .id = "T") %>% 
  mutate(D = ifelse(briefwahlbezirk %in% treated & substr(briefwahlbezirk,1,2) == "12" ,1,0),
         `T` = (as.numeric(`T`) - 2)*-1) 

stimmen <- c("Erststimme", "Zweitstimme")

plot_data <- bind_rows(map(stimmen, get_did_counterfactual))

did_data <- data_2024 %>% 
  select(id, partei, stimme, sim) %>% 
  rename(briefwahlbezirk = id,
         anteil = sim) %>% 
bind_rows(select(data, briefwahlbezirk, partei, stimme, anteil),
          .id = "T") %>% 
  mutate(D = ifelse(briefwahlbezirk %in% treated & substr(briefwahlbezirk,1,2) == "12" ,1,0),
         `T` = (as.numeric(`T`) - 2)*-1) 

stimmen <- c("Erststimme", "Zweitstimme")

plot_data <- bind_rows(map(stimmen, get_did_counterfactual))

plot_data %>% 
  ggplot(aes(`T`,Y,color = D)) +
  scale_y_continuous("Stimmenanteil", labels = scales::label_percent(1)) + 
  scale_x_continuous("Wahlen", limits = c(20.8,25), breaks = c(21,24), labels = function(x) paste0("20",x)) +
  scale_color_manual(NULL, values = c(paste0("#0E8C1D","88"), "#0E8C1D", "darkgrey")) +
  scale_linetype_manual(NULL, values = c("solid","solid","dashed")) +
  geom_line(aes(group = D, linetype = D)) +
  geom_point(size = 1) +
  ggbrace::stat_brace(data = bind_rows(map(stimmen, get_brace_data)), aes(x,y), rotate = 90, color = "black", outside = F) +
  geom_text(data = bind_rows(map(stimmen, get_text_data)), aes(x,y, label = lab), hjust = -0.1, color = "black", size = 4) +
  facet_wrap(~stimme) +
  theme_light(base_size = 18) +
  theme(legend.position = "bottom",
        panel.grid = element_blank()) #+
  ggtitle('Difference-in-Differences-Schätzung des "Plakateffekts"')

ggsave(file = "did_erst.png", width = 25, height = 15, units = "cm", dpi = 800)


fit_feols <- function(p, s) feols(anteil ~ D*`T`, data = filter(did_data, partei == p & stimme == s), cluster = "briefwahlbezirk") 

parteien <- c("SPD", "CDU", "GRÜNE", "AfD", "FDP", "DIE LINKE", "Sonstige")

map2(c(parteien,parteien), c(rep(stimmen[1],7),rep(stimmen[2],7)), fit_feols) %>% 
  map(~broom::tidy(.x, conf.int = T)) %>% 
  bind_rows(.id = "ps") %>% 
  mutate(ps = as.numeric(ps),
         partei = factor(parteien[ifelse(ps > 7,ps - 7 ,ps)],
                         levels = rev(parteien)),
         stimme = stimmen[ifelse(ps > 7,2,1)]) %>% 
  filter(term == "D:T") %>% 
  ggplot(aes(estimate, partei,xmin = conf.low, xmax = conf.high, color = partei)) +
  scale_x_continuous("ATT Schätzungen", labels = function(x) paste0(round(x*100,0),"%P")) +
  scale_y_discrete(NULL) +
  scale_color_manual(values = rev(c("#CC0000", "#000000", "#0E8C1D", "#005EA4", "#FFC000", "#CC0066", "#BFBFBF")), guide = "none") +
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_pointrange(size = 1, linewidth = 1) +
  facet_wrap(~stimme) +
  theme_light(base_size = 18) +
  theme(panel.grid = element_blank())

ggsave(file = "did_complete.png", width = 25, height = 15, units = "cm", dpi = 800)

# Add ons ----

agh_geo <- st_read(here("data","RBS_OD_UWB_AH21", "RBS_OD_UWB_AH21.shp")) %>% 
  rename(postal_id = BWB) %>% 
  mutate(wdh = ifelse(postal_id %in% treated, "Wiederholung", "Keine Wiederholung"),
         rdf = factor(ifelse(BWK == "77", 1, 0)))

bernstein_colors <-  c("#003255", "#A02D23", "#378981", "#415F7D", "#A5B4B4", "#601B15", "#6E7887", "#21524D", "#464B4B", "#27394B")

ggplot(agh_geo, aes(fill = wdh, geometry = geometry, alpha = rdf)) +
  scale_fill_manual(NULL, values = c("grey", bernstein_colors[2])) +
  scale_alpha_manual(values = c(0.6, 1), guide = "none") +
  geom_sf(color = "white", linewidth = 0.05) +
  ggtitle("Wiederholung BTW im Wahlkreis 77") +
  theme_void() +
  theme(legend.position = "bottom") 

ggsave(file = "map.png", width = 16, height = 9, units = "cm", dpi = 600)

agh_geo %>% 
  left_join(data_matched, by = c("postal_id" = "id")) %>% 
  mutate(D = factor(case_when(is.na(D) & wdh == "Keine Wiederholung" ~ wdh,
                              is.na(D) & wdh != "Keine Wiederholung" ~ "Ausgeschlossen",  
                              BWK == "77" & !is.na(D) ~ "Untersuchungsgruppe",
                              BWK != "77" & !is.na(D) ~ "Kontrollgruppe"))) %>% 
  ggplot(aes(fill = D, geometry = geometry, alpha = wdh)) +
  scale_fill_manual(NULL, values = c("grey", "grey", bernstein_colors[1], bernstein_colors[2])) +
  scale_alpha_manual(values = c(0.6, 1), guide = "none") +
  geom_sf(color = "white", linewidth = 0.05) +
  guides(fill=guide_legend(ncol = 2)) +
  theme_void() +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white", color = "white")) 

ggsave(file = "map_tc.png", width = 16, height = 9, units = "cm", dpi = 600)

social %>% 
  mutate(weight = ifelse(id %in% data_matched$id,1,0)) %>% 
  tidyr::pivot_longer(-c(id,D,weight)) %>%
  mutate(value = value/100) %>% 
  calc_coefs() %>% 
  filter(!name %in% c("AG", "Deutsche18VeränderungzumVo")) %>% 
  rename(Variables = name) %>%  
  tidyr::pivot_longer(-c(Variables,df)) %>% 
  mutate(Matched = ifelse(str_detect(name, "w_"),"Mit Matching", "Ohne Matching"),
         Statistic = ifelse(str_detect(name, "est"), "est", "se")) %>% 
  dplyr::select(-name) %>% 
  tidyr::pivot_wider(names_from = "Statistic", values_from = "value") %>% 
  mutate(t_crit = qt(0.025, df, lower.tail = F),
         lwr = est - t_crit*se,
         upr = est + t_crit*se,
         in_eb = ifelse(Variables %in% colnames(social)[bal_vars],"Ja","Nein"),
         Gs = factor(case_when(str_detect(Variables, "2021") ~ "Wahlverhalten",
                               Variables %in% c("Einwohnerunter6Jahre", "Einwohner6bisunter18", "Einwohner18bisunter65", "Einwohner65") ~ "Alter",
                               Variables %in% c("Deutsche18bisunter25", "Deutsche25bisunter35", "Deutsche35bisunter45", "Deutsche45bisunter60", "Deutsche60bisunter70", "Deutsche70") ~ "Wahlberechtigte",
                               Variables %in% c("Deutsche18Familienstandledi", "Deutsche18Familienstandverh", "Deutsche18Familienstandverw", "Deutsche18Familienstandgesc", "Deutsche18FamilienstandELP") ~ "Familienstand",
                               Variables %in% c("Deutsche16mitMigrationshint", "Deutsche18mitMigrationshint", "Einwohnerunter65inSGBII201", "Deutsche18weiblich", "Deutsche18Religionszugehörig", "Ausländer", "EUBürger16", "Deutsche16", "Deutsche18") ~ "Andere"),
                     levels = c("Wahlverhalten", "Wahlberechtigte", "Alter", "Familienstand", "Andere")),
         Variables = case_when(Gs == "Wahlverhalten" ~ str_remove(Variables,"2021"),
                               Gs == "Wahlberechtigte" ~ str_replace_all(str_replace_all(str_replace_all(Variables,"Deutsche", "Deutsche "),"bis", " bis "), "unter", "unter "),
                               Gs == "Alter" ~ str_replace_all(str_replace_all(str_replace_all(str_replace(Variables, "Jahre", " Jahre"),"Einwohner", "Einwohner "),"bis", " bis "), "unter", "unter "),
                               Gs == "Familienstand" ~ str_remove_all(Variables, "Deutsche18Familienstand"),
                               T ~ Variables)) %>% 
  filter(Gs != "Andere") %>% 
  ggplot(aes(est, Variables, color = Matched)) +
  scale_y_discrete(NULL) +
  scale_x_continuous("Unterschied zwischen Untersuchungs- und Kontrolgruppe", labels = function(x) paste0(round(x*100),"%P")) +
  scale_color_manual(NULL, values = c(bernstein_colors[3], bernstein_colors[4])) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(xmin = lwr, xmax = upr), position = position_dodge(0.5), alpha = 0.2) +
  guides(color=guide_legend(ncol = 2)) +
  facet_grid(rows = vars(Gs), scales = "free_y") +
  theme_light(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

ggsave(file = "tc.png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave(file = "file_name.png", width = 31, height = 12, units = "cm", dpi = 600)

wdh_data %>% 
  filter(partei == "GRÜNE" & stimme == "Zweitstimme")  %>% 
  left_join(select(data_2024, id, partei, gueltige_stimmen, stimme), by = c("briefwahlbezirk" = "id", "partei", "stimme")) %>% 
  mutate(abs_21 = anteil * gultige_stimmen,
         abs_24 = sim * gueltige_stimmen) %>% 
  summarise(across(c(gultige_stimmen, gueltige_stimmen, abs_21, abs_24), \(x) sum(x)), .by = D) %>% 
  mutate(p_21 = abs_21/gultige_stimmen,
         p_24 = abs_24/gueltige_stimmen,
         D = ifelse(D == 1, "Reinickendorf", "Rest Berlin")) %>% 
  select(D, abs_21:p_24) %>% 
  tidyr::pivot_longer(cols = -D, names_pattern = "(.*)_(.*)", names_to = c("type", "year")) %>% 
  tidyr::pivot_wider(names_from = type, values_from = value) %>% 
  ggplot(aes(year, p)) +
  scale_y_continuous("Zweitstimmenanteil", labels = scales::label_percent(1), expand = c(0.001, 0.001), limits = c(0, 0.30)) +
  scale_x_discrete(NULL, labels = function(x) paste0("20",x)) +
  geom_col(fill = "#0E8C1D") +
  geom_text(aes(label = paste0(round(p*100,1),"%")), vjust = 1.5, color = "white", hjust = 0.5, size = 5) +
  theme_light(base_size = 18) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~D)


ggsave(file = "baseline_diffs.png", width = 33.9, height = 19.1, units = "cm", dpi = 800)



