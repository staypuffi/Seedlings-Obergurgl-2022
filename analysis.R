### Obergurgl seedlings Analysis ----

# run esploration script to get sum tables
source("exploration.R")

library(lme4)
library(lmerTest)  # mixed models
library(sjPlot)  # to get model tables
library(vegan)
library(modEvA)  # for variation partiting

# Regress overall germination rates (=total number of seedlings) on the 
 # microclimate of the experimental plots (suggested: summer temperatures (June/July) , 
 # alternatively also: length of snow cover – simple Poisson regression with numbers or 
 # binomial one with rates.

### total germination rates over Temperature & snow_cover (GLM) ----
sum_temp
all_seeds

## over all species summed together
  area_seeds <- all_seeds %>%
    left_join(sum_temp) %>%
    group_by(site, area, mowing) %>%
    summarise(tot_seedlings = sum(number_seedlings),
              germ_rate = tot_seedlings/(20*50) ) %>%
    left_join(sum_temp)
  
  area_seeds_persp <- seeds_and_sumtemp %>%
    group_by(site, area, mowing, taxon) %>%
    summarise(seedlings = sum(number_seedlings),
              germ_rate = seedlings/20 ) %>%
    left_join(sum_temp)
  
  ## plot 
    # temperature
    area_seeds %>%
      ggplot(aes(x= summer_mean, y= tot_seedlings,
                 col= mowing))+
      geom_point()+
      theme_classic()+
      geom_smooth(method = "glm",
                  method.args = list(family = "poisson"))
    # snow
    area_seeds %>%
      ggplot(aes(x= snow_days, y= tot_seedlings,
                 col= mowing))+
      geom_point()+
      theme_classic()+
      geom_smooth(method = "glm",
                  method.args = list(family = "poisson"))
  
    # snowmelt
    area_seeds %>%
      ggplot(aes(x= last_snowday, y= tot_seedlings,
                 col= mowing))+
      geom_point()+
      theme_classic()+
      geom_smooth(method = "glm",
                  method.args = list(family = "poisson"))
    
    
    ## test
    area_seeds <- area_seeds %>%
      mutate(summer_mean_from0= summer_mean - min(summer_mean), # start model from 0
              snow_days_from0= snow_days - min(snow_days),
             last_snowday_from0= last_snowday - min(last_snowday))  
    
    germ_temp_mod <- glm(data = area_seeds,
                tot_seedlings ~ summer_mean_from0 +
                      snow_days_from0 +
                      last_snowday_from0+
                      mowing,
                family = "poisson")
    
    summary(aov(germ_temp_mod))
    summary(germ_temp_mod)
  
    plot(germ_temp_mod)

  ## plot model 
  model_pred <- as.data.frame(predict(germ_temp_mod)) 
  model_pred <- model_pred %>%
    mutate(plot_ID = rownames(model_pred)) %>%
    rename("predictions" = "predict(germ_temp_mod)") %>%
    as_tibble()
  
 ## variation partitioning
  # calulate the variance explained by the complete model,  
  # then calculate models with all variables seperateky and also D², =the variance explained
  # then try to partition following a Venn diagram
 
  Tot <- with(summary(germ_temp_mod), 1 - deviance/null.deviance)
 
 # Temp = A, Snowtime = B, Snowmelt = C, site = D
 
 ## 3 vars
 ABC <- with(summary(glm(data = area_seeds,
                               tot_seedlings ~ summer_mean_from0+
                                                snow_days_from0+
                                                last_snowday_from0,
                               family = "poisson")), 1 - deviance/null.deviance)
 ABD <- with(summary(glm(data = area_seeds,
                                 tot_seedlings ~ summer_mean_from0+
                                 snow_days_from0+
                                 site,
                                family = "poisson")), 1 - deviance/null.deviance)
 ACD <- with(summary(glm(data = area_seeds,
                              tot_seedlings ~ summer_mean_from0+
                                last_snowday_from0 +
                                site,
                              family = "poisson")), 1 - deviance/null.deviance)
 BCD <- with(summary(glm(data = area_seeds,
                              tot_seedlings ~ snow_days_from0+
                                last_snowday_from0 +
                                site,
                              family = "poisson")), 1 - deviance/null.deviance)
 ## 2 vars
 AB <- with(summary(glm(data = area_seeds,
                              tot_seedlings ~ summer_mean_from0+
                                snow_days_from0,
                              family = "poisson")), 
                 1 - deviance/null.deviance)
 AC <- with(summary(glm(data = area_seeds,
                             tot_seedlings ~ summer_mean_from0+
                               last_snowday_from0,
                             family = "poisson")), 1 - deviance/null.deviance)
 BC <- with(summary(glm(data = area_seeds,
                        tot_seedlings ~ snow_days_from0+
                          last_snowday_from0,
                        family = "poisson")), 1 - deviance/null.deviance)
 AD <- with(summary(glm(data = area_seeds,
                        tot_seedlings ~ summer_mean_from0+
                          site,
                        family = "poisson")), 1 - deviance/null.deviance)
 BD <- with(summary(glm(data = area_seeds,
                             tot_seedlings ~ snow_days_from0+
                               site,
                             family = "poisson")), 1 - deviance/null.deviance)
 CD <- with(summary(glm(data = area_seeds,
                             tot_seedlings ~ last_snowday_from0+
                               site,
                             family = "poisson")), 1 - deviance/null.deviance)
 ## only 1 var
 A <- with(summary(glm(data = area_seeds,
     tot_seedlings ~ summer_mean_from0 ,
     family = "poisson")), 1 - deviance/null.deviance)
 B <- with(summary(glm(data = area_seeds,
        tot_seedlings ~ snow_days_from0 ,
        family = "poisson")), 1 - deviance/null.deviance)
 C <- with(summary(glm(data = area_seeds,
              tot_seedlings ~ last_snowday_from0 ,
              family = "poisson")), 1 - deviance/null.deviance)
 D <- with(summary(glm(data = area_seeds,
                tot_seedlings ~ site ,
                 family = "poisson")), 1 - deviance/null.deviance)
 ## Venn DIagram
 Tot
 ## 2 vars
 #AB
 AoB <- AB-B
 BoA <- AB-A
 AuB <- AB-(AoB+BoA)
 #AC
 AoC <- AC-C
 CoA <- AC-A
 AuC <- AC-(AoC+CoA)
 #AD
 AoD <- AD-D
 DoA <- AD-A
 AuD <- AD-(AoD+DoA)
 #BC
 BoC <- BC-C
 CoB <- BC-B
 BuC <- BC-(BoC+CoB)
 #BD
 BoD <- BD-D
 DoB <- BD-B
 BuD <- BD-(BoD+DoB)
 #CD
 CoD <- CD-D
 DoC <- CD-C
 CuD <- CD-(CoD+DoC)
 ## 3 vars
 
 AoBoC <- ABC-B-C+BuC
 AoBoD <- ABD-C-D+CuD
 AoCoD <- ACD-C-D+CuD
 
 BoAoC <- ABC-A-C+AuC
 BoCoD <- ABD-C-D+CuD
 BoAoD <- ABD-A-D+AuD
 
 CoAoB <- ABC-A-B+AuB
 CoBoD <- ABD-B-D-BuD
 CoAoD <- ACD-A-D+AuD
 
 DoAoB <- ABD-A-B+AuD
 DoAoC <- ACD-A-C+AuD
 DoBoC <- BCD-B-C+BuC
 
 AuBuC <- ABC-A-B-C+AuB+AuC+BuC
 AuBuD <- ABD-A-B-D+AuB+AuD+BuD 
 AuCuD <- ACD-A-C-D+AuC+AuD+CuD
 BuCuD <- BCD-B-C-D+BuC+BuD+CuD
 
 ## 4 vars
 AuBuCuD <- Tot-A-B-C-D+AuB+AuC+AuD+BuC+BuD+CuD-AuBuC-AuBuD-AuCuD-BuCuD
 
 # lei A
 Tot-B-C-D+AuB+AuC+AuD+BuC+BuD+CuD-AuBuC-AuBuD-AuCuD-BuCuD-AuBuCuD
 ### oder eaaasy lösung...
 (Tot-BCD)/Tot*100
 (Tot-ACD)/Tot*100
 (Tot-ABD)/Tot*100
 (Tot-ABC)/Tot*100
 
 var_expl <- data.frame(variable = c("overall", "summer_temp", "snow_duration", "snowmelt_date", "mowing"),
                        var_explained=c(Tot, Tot-BCD, Tot-ACD, Tot-ABD, Tot-ABC))
 var_expl <- var_expl %>%
   mutate(var_expl_prop = round(var_explained/Tot*100,1),
          var_explained = round(var_explained, 4))
 
 
 ggvenn::ggvenn(show_elements = T,
   data = list('Summer temp'=round(c(var_expl[2,2], AuB-AuBuC-AuBuCuD, AuC-AuCuD-AuBuCuD, AuD-AuCuD-AuBuCuD, AuCuD-AuBuCuD, AuBuC-AuBuCuD, AuBuD-AuBuCuD, AuBuCuD)/Tot, 3),
               'Snow duration'=round(c(var_expl[3,2], AuB-AuBuC-AuBuCuD, BuC-BuCuD-AuBuCuD, BuD-BuCuD-AuBuCuD, AuBuC-AuBuCuD, AuBuD-AuBuCuD, BuCuD-AuBuCuD, AuBuCuD)/Tot, 3),
              'Snowmelt date'=round(c(var_expl[4,2], AuC-AuCuD-AuBuCuD, BuC-BuCuD-AuBuCuD, CuD-AuCuD-AuBuCuD, BuCuD-AuBuCuD, AuBuC-AuBuCuD, AuCuD-AuBuCuD, AuBuCuD)/Tot, 3),
               'Mowing'=round(c(var_expl[5,2], AuD-AuCuD-AuBuCuD, BuD-BuCuD-AuBuCuD, CuD-AuCuD-AuBuCuD, BuCuD-AuBuCuD, AuCuD-AuBuCuD, AuBuD-AuBuCuD, AuBuCuD)/Tot, 3))
 )

### derive natural frequency of species at sites to use as co variable (to account for bias) ----
      unique(gurgl_species$cover)
      
      all_seeds <- gurgl_species %>%
        mutate(cover_perc = case_when(cover == "+" |
                                        cover == "p" | cover == "r" ~ 1,
                                      cover == "1a" ~ 3,
                                      cover == "1b" ~ 5,
                                      cover == "2a" ~ 10,
                                      cover == "2b" ~ 15,
                                      cover == "3" ~ 25,
                                      cover == "4" ~ 40,
                                      cover == "5" ~ 70,
                                      TRUE ~ as.double(cover) )) %>% 
        select(taxon:mowing, cover, cover_perc) %>%
        filter(!is.na(taxon)) %>%
        group_by(taxon, site) %>% # now reduce to mean frequency on site
        summarise(natural_freq = mean(cover_perc)) %>%
        right_join(all_seeds) %>% # join with seed data
        mutate(natural_freq = replace_na(natural_freq, 0)) %>%
        relocate(natural_freq, .after = germ_rate)
      
      
### derive species niches -----
schran_niches <- schran_species %>%
  left_join(schran_sumtemp)

(schran_niches <- schran_niches %>%
  group_by(TAXON) %>%
    summarise(summer_niche = median(summer_mean, na.rm = T),
              snow_niche = median(snow_days, na.rm = T),
              snowmelt_niche = median(last_snowday, na.rm = T)) %>%
  rename("taxon" = "TAXON")
)

# combine niches with seedling number
      write_csv2(all_seeds %>%
        group_by(taxon) %>%
        summarise(number_plots_sown = length(taxon)) %>%
        right_join(tot_seeds_species) %>%
        mutate(germ_rate = (number_seedlings/(number_plots_sown*20))*100) %>%
        left_join(schran_niches) %>%
          select(taxon, summer_niche, number_seedlings, number_plots_sown, germ_rate),
        "plots/species_niches_and_germrate.csv")
      
      
      
### niche diff for each plot (mixed model) ----
  
niche_seeds <- all_seeds %>%
  left_join(schran_niches) %>%
  left_join(sum_temp) %>%
  mutate(dist_summer = summer_mean-summer_niche,
         dist_snow = snow_days-snow_niche,
         dist_snowmelt = last_snowday-snowmelt_niche) %>%
  relocate(summer_mean, .after = summer_niche) %>%
      relocate(dist_summer, .after = summer_mean) %>%
  arrange(site, area, taxon)

      niche_seeds %>% select(site, area, summer_niche, 
                             summer_mean, dist_summer)

## plotting niches
  # with original data to see original spread
  schran_species %>%
    filter(TAXON %in% my_spec$taxon) %>%
     left_join(schran_sumtemp) %>%
      rename("taxon" = "TAXON") %>%
      left_join(schran_niches) %>%
        ggplot(aes(x= reorder(taxon, snow_niche),
                   y= snow_days))+
        geom_boxplot()+
            #ylim(c(4,17))+
        coord_flip()+
        theme_minimal()      
      
niche_seeds %>%
  ggplot(aes(x= reorder(taxon, snowmelt_niche),
             y= snowmelt_niche))+
  geom_point()+
  coord_flip()+
  theme_minimal()

niche_seeds %>%
  ggplot(aes(x= summer_niche,
             y= snowmelt_niche))+
  geom_point()+
  theme_minimal()
  
  # distance from niche at site and germ rate

  niche_seeds %>%
    ggplot(aes(x= dist_summer, y= number_seedlings,
               col = mowing
               ))+
    geom_point()+
    geom_smooth(method = "glm",
                formula = "y ~ poly(x, 2)",
                method.args = list(family = "poisson"),
                col = "black", linetype = 2, size = .75
                )+
   geom_vline(xintercept = 0, linetype = 2)+
    facet_wrap(~site, scales = "free_x")+
    theme_classic()+
    scale_shape_manual(values = c(1,16))

  niche_seeds %>%
    ggplot(aes(x= natural_freq, y= number_seedlings))+
    geom_point()+
  geom_smooth(method = "glm")
    

## tests
  min(niche_seeds$dist_summer)

  niche_distoptimum_glmer <- glmer(data = niche_seeds, 
      number_seedlings ~ dist_summer + mowing +
                                  natural_freq +
                                  (1|site/area),
            family = "poisson")

exp(coef(niche_distoptimum_glmer)[[1]])
summary(niche_distoptimum_glmer)

### variation partitoning
  #res <- partR2::partR2(data = niche_seeds,
   #     mod =niche_distoptimum_glmer,
    #    partvars = c("dist_summer", "mowing", 
     #                #"dist_summer:mowing", 
      #               "natural_freq"), nboot=10)

                      
#### regress niches over overall germ rate (6 GLMs) ----

    # start temperatures at min-value
    min_summerniche <- min(niche_seeds$summer_niche)
    
niche_seeds <- niche_seeds %>%
      mutate(summer_niche_from0 = summer_niche-min_summerniche) 

## combine to site level to avoid replication

niche_seeds_overall <- niche_seeds %>%
  group_by(site, taxon) %>%
  summarise(summer_niche = mean(summer_niche, na.rm = T),
            summer_niche_from0 = mean(summer_niche_from0, na.rm = T),
            number_seedlings = sum(number_seedlings, na.rm = T),
            germ_rate = number_seedlings/(20*5),
            natural_freq = mean(natural_freq, na.rm = T))

niche_seeds_allsites <- niche_seeds_overall %>%
  group_by(taxon) %>%
  summarise(summer_niche = mean(summer_niche, na.rm = T),
            summer_niche_from0 = mean(summer_niche_from0, na.rm = T),
            number_seedlings = sum(number_seedlings, na.rm = T),
            germ_rate = number_seedlings/(20*5),
            natural_freq = mean(natural_freq, na.rm = T))


## plot
niche_seeds_overall %>%
  ggplot(aes(x= summer_niche, y = number_seedlings))+
  geom_point()+
  #geom_text(aes(label = taxon))+
  geom_smooth(method = "glm",
              method.args = list(family = "poisson"),
              formula = "y ~ x")+
  theme_classic()+
  facet_wrap(~site)

niche_seeds_allsites %>%
  ggplot(aes(x= summer_niche, y = number_seedlings))+
  geom_point()+
  #geom_text(aes(label = taxon))+
  geom_smooth(method = "glm",
              method.args = list(family = "poisson"),
              formula = "y ~ x")

## test
model_coefs <- tibble(site = c("all", 1:5),
                      Intercept = NA, 
                      slope = NA)
  # all sites
  niche_germrate_all <- glm(
        data = niche_seeds_allsites,
          number_seedlings ~ summer_niche_from0 + natural_freq,
        family = "poisson")
  
  summary(niche_germrate_all)
  summary(aov(niche_germrate_all))
  model_coefs[1, 2] <- coef(niche_germrate_all)[1]
  model_coefs[1, 3] <- coef(niche_germrate_all)[1]+coef(niche_germrate_all)[2]
  
  #site 1
  niche_seeds_sitespecif <- niche_seeds_overall %>%
    filter(site == "1")
  
  niche_germrate_1 <- glm(
    data = niche_seeds_sitespecif,
    number_seedlings ~ summer_niche_from0 + natural_freq,
    family = "poisson")
  
  summary(niche_germrate_1)
  summary(aov(niche_germrate_1))
  model_coefs[2, 2] <- coef(niche_germrate_1)[1]
  model_coefs[2, 3] <- coef(niche_germrate_1)[1]+coef(niche_germrate_1)[2]
  
  
   #site 2
  niche_seeds_sitespecif <- niche_seeds_overall %>%
    filter(site == "2")
  
  niche_germrate_2 <- glm(
    data = niche_seeds_sitespecif,
    number_seedlings ~ summer_niche_from0 + natural_freq,
    family = "poisson")
  
  summary(niche_germrate_2)
  summary(aov(niche_germrate_2))
  model_coefs[3, 2] <- coef(niche_germrate_2)[1]
  model_coefs[3, 3] <- coef(niche_germrate_2)[1]+coef(niche_germrate_2)[2]
  
   #site 3
  niche_seeds_sitespecif <- niche_seeds_overall %>%
    filter(site == "3")
  
  niche_germrate_3 <- glm(
    data = niche_seeds_sitespecif,
    number_seedlings ~ summer_niche_from0 + natural_freq,
    family = "poisson")
  
  summary(niche_germrate_3)
  summary(aov(niche_germrate_3))
  model_coefs[4, 2] <- coef(niche_germrate_3)[1]
  model_coefs[4, 3] <- coef(niche_germrate_3)[1]+coef(niche_germrate_3)[2]
  
  #site 4
  niche_seeds_sitespecif <- niche_seeds_overall %>%
    filter(site == "4")
  
  niche_germrate_4 <- glm(
    data = niche_seeds_sitespecif,
    number_seedlings ~ summer_niche_from0 + natural_freq,
    family = "poisson")
  
  summary(niche_germrate_4)
  summary(aov(niche_germrate_4))
  model_coefs[5, 2] <- coef(niche_germrate_4)[1]
  model_coefs[5, 3] <- coef(niche_germrate_4)[1]+coef(niche_germrate_4)[2]
  
  #site 5
  niche_seeds_sitespecif <- niche_seeds_overall %>%
    filter(site == "5")
  
  niche_germrate_5 <- glm(
    data = niche_seeds_sitespecif,
    number_seedlings ~ summer_niche_from0 + natural_freq,
    family = "poisson")
  
  summary(niche_germrate_5)
  summary(aov(niche_germrate_5))
  model_coefs[6, 2] <- coef(niche_germrate_5)[1]
  model_coefs[6, 3] <- coef(niche_germrate_5)[1]+coef(niche_germrate_5)[2]
  
  ## summary
  model_coefs
  cbind(model_coefs[,1], round(exp(model_coefs[,2:3]),2))
  
  summary(aov(niche_germrate_all))
  summary(aov(niche_germrate_1))
  summary(aov(niche_germrate_2))
  summary(aov(niche_germrate_3))
  summary(aov(niche_germrate_4))
  summary(aov(niche_germrate_5))

  ## Variation partitioning
  varPart(model = niche_germrate_all,
          groups = as.data.frame(cbind(
            c("summer_niche_from0", "natural_freq"),
            rep("niche", 2) ))
          )

#### extra: germination niches ------
  all_seeds %>% filter (number_seedlings > 1) %>% select(taxon) %>% unique()
  
seed_niches <- all_seeds %>% ungroup() %>%
    left_join(sum_temp) %>%
    filter(number_seedlings > 1) %>%
    group_by(taxon) %>%
    summarise(seed_summer_niche = median(summer_mean, na.rm = T),
              seed_snow_niche = median(snow_days, na.rm = T),
              seed_summer_rangemin = min(summer_mean, na.rm = T),
              seed_summer_rangemax = max(summer_mean, na.rm = T)) %>%
    left_join(schran_niches)
  
## overview of seed niches

  all_seeds %>%
    left_join(sum_temp) %>%
    filter(number_seedlings > 1) %>%
    ggplot(aes(x=snow_days, y= number_seedlings,
               col = taxon))+
    geom_point(show.legend = FALSE)+
    facet_wrap(~taxon, scales = "free_y")+
    theme_classic()+
    geom_smooth(method = "glm",
                method.args = list(family = "gaussian"),
                formula = "y ~ poly(x, 2)",
                se = F,
                show.legend = FALSE)
  
# do they correlate with adult niches from Schrankogel?
  seed_niches %>%
    ggplot(aes(x= snow_days, y=seed_snow_niche))+
    geom_point()+
    geom_smooth(method = "glm")
      
# with number seedlings (= germ success)?
 tot_seeds_species <- seed_niches %>%
    right_join(tot_seeds_species) 
  
 tot_seeds_species %>%
   ggplot(aes(x= seed_snow_niche, y=number_seedlings))+
    geom_point()+
    geom_smooth(method = "glm",
                method.args = list(family = "poisson"))
  
   
    summary(glm(data = tot_seeds_species, 
        number_seedlings ~ seed_snow_niche,
        family = "poisson"))
  
    