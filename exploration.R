##### Obergurgl seeds #######
## 5.9.2022
library(tidyverse)
library(ggplot2)
library(lubridate)
setwd("~/Uni/Mec 8 - Obergurgl seeds/Obergurgl seeds")

#### load data ####
# seeds
all_seeds <- read_csv("data/seedlings2022_obergurgl_correct.csv")
all_seeds

(my_spec <- unique(all_seeds$taxon))
write_csv2(as_tibble(my_spec), "data/seedling_species.csv")
my_spec <- read_csv2("data/seedling_species_chars.csv")

table(my_spec$growthform)


# species
schran_species <- read_csv2("data/Schrankogel_species_data.csv")
schran_species

gurgl_species <- read_csv2("data/Obergurgl_species_data.csv")
gurgl_species


## check wich species s´missing in seedlings
  # schran (for niches)
    setdiff(all_seeds$taxon, schran_species$TAXON)
    sort(setdiff(schran_species$TAXON, (all_seeds$taxon)))
  
  # gurgl (for bias control)
    setdiff(all_seeds$taxon, gurgl_species$taxon)
    sort(setdiff(gurgl_species$taxon, all_seeds$taxon))
    
# Desch cesp and Carex curv have subsp. in SChran but not in Obergurgl, remove for compatibility
schran_species <- schran_species %>%
  mutate(TAXON = case_when(TAXON == "Carex curvula subsp. curvula" ~ "Carex curvula",
                           TAXON == "Deschampsia cespitosa subsp. cespitosa" ~ "Deschampsia cespitosa",
                           TRUE ~ TAXON))
gurgl_species <- gurgl_species %>%
  mutate(taxon = case_when(taxon == "Carex curvula subsp. curvula" ~ "Carex curvula",
                           taxon == "Deschampsia cespitosa subsp. cespitosa" ~ "Deschampsia cespitosa",
                           TRUE ~ taxon))



#### Temp loggers ----
gurgl_temp <- read_csv("data/2022-11-30_Loggers_readout_Obergurgl.csv")
gurgl_temp

unique(gurgl_temp$block)

gurgl_temp <- gurgl_temp %>%
              filter(block == "A") %>%
              rename(mowing = treatment)%>%
              arrange(site, area, mowing, date) %>%
              mutate(day = date(date),
                     month = month(date),
                     hour = hour(date),
                     site = as.factor(site),
                     area = as.factor(area),
                     mowing = as.factor(mowing))



levels(gurgl_temp$site) ## site 2 completely missing
levels(gurgl_temp$area) ## site 2 completely missing
levels(gurgl_temp$mowing)

## Schrankogel temps -
schran_temp <- read_csv("data/2022-11-27_Loggers_readout_Schrankogel.csv")
schran_temp

schran_temp <- schran_temp %>%
  mutate(day = date(date),
         month = month(date),
         hour = hour(date))

 ### split Obergurgl plot column into Site, Area, Block and mowing -----
  all_seeds %>%
    filter(!is.na(number_seedlings)) %>%
    separate(plot, c("site", "area", "block", "mowing"),
             sep = "-") %>%
    mutate(germ_rate = number_seedlings/20) %>%
    select(-block) -> all_seeds # we can delete blocks as we only have block A

  gurgl_species <- gurgl_species %>%
    separate(plot_ID, c("site", "area", "block", "mowing"),
             sep = "-") %>%
    mutate(across(.cols = 1:4, .fns=factor),
           taxon = factor(taxon, levels = unique(gurgl_species$taxon))) 
  
  
  ## which species are not present in all plots? (unfinished)
  table(all_seeds$taxon)[table(all_seeds$taxon) != 50]
    
    #so we have to exclude Calluna vulgaris from the analysis (to make the plots comparable)  
  all_seeds <- all_seeds %>% filter(taxon != "Calluna vulgaris")
  
  table(all_seeds$taxon)
  
## total seeds
  tot_seeds <- all_seeds %>%
    group_by(site, area, mowing) %>%
    summarise(tot_seeds = sum(number_seedlings),
              germ_rate = sum(germ_rate))
# and per spec
  all_seeds %>% 
    group_by(taxon) %>% 
    summarise(number_seedlings = sum(number_seedlings),
              germ_rate = sum(germ_rate)) %>%
    arrange(-number_seedlings) -> tot_seeds_species 
  
  
#### summarize temps ----
  # Obergurgl
    # do we have all data?
  gurgl_temp %>%
    group_by(site, area) %>%
    summarise(what_mowing_is_there = paste(unique(mowing)))
  
# annual
  sum_temp <- gurgl_temp %>%
    group_by(site, area, mowing) %>%
    summarise(ann_mean = mean(temp))

# summer
  sum_temp <- gurgl_temp %>%
    filter(month == 6 | month == 7) %>%
    group_by(site, area, mowing) %>%
    summarise(summer_mean = mean(temp)) %>%
    full_join(sum_temp)

# number of snow days
  ## test
  sum_temp <- gurgl_temp %>%
    group_by(site, area, mowing, day) %>%
    summarise(temp_mean = mean(temp),
              temp_min = min(temp),
              temp_max = max(temp)) %>%
    filter(#temp_mean < 1 & temp_mean > -1 &
             temp_min > -1 &
             temp_max < 1 ) %>%
    #mutate(month = month(day)) %>%
    group_by(site, area, mowing ) %>% #, month) %>%
      summarise(snow_days = length(temp_mean)) %>%
    full_join(sum_temp) %>%
    mutate(snow_days = replace_na(snow_days, 0) )# one site has NA snow
  
  
  ## plot with month to get a feeling how good the calculation is
  sum_temp %>%
    ggplot(aes(y=snow_days, x=site, fill = site))+
    geom_boxplot()
  
  sum_temp %>%
    ggplot(aes(y=snow_days, x=month, fill = site))+
    geom_col()+
  scale_x_continuous(breaks = seq(1,12, by = 1))
  
  sum_temp %>%
    ggplot(aes(y=snow_days, x=month, col = area, 
               pch = area))+
    geom_line(type = 1)+
    facet_wrap(~site)+
    theme_classic()
    
      # max min between +/- 1 optimal
  
  (sum_temp <- sum_temp %>% arrange(site, area))
  
  ## test different ways of calc snowcover
  sum_temp %>%
  ggplot(aes(y=snow_days, x=summer_mean,
             col= site, pch = mowing, linetype = mowing))+
    geom_point(aes())+
    geom_smooth(se=F,
                method = "glm")
  
 ## day of snowmelt
  sum_temp <- gurgl_temp %>%
    group_by(site, area, mowing, day) %>%
    summarise(temp_mean = mean(temp),
              temp_min = min(temp),
              temp_max = max(temp)) %>%
    filter(#temp_mean < 1 & temp_mean > -1 &
      temp_min > -1 &
        temp_max < 1) %>%
    group_by(site, area, mowing) %>%
      summarise(last_snowday = max(day)) %>%
    right_join(sum_temp)
  
   sum_temp %>%
     ggplot(aes(x= site, y= last_snowday))+
     geom_boxplot()
  
   
# number of frost days per month
   frost_days <- gurgl_temp %>%
     group_by(site, area, mowing, day) %>%
     summarise(temp_min = min(temp)) %>%
     filter(temp_min < -0.5 &
              day > "2022-04-30") %>%
     mutate(month = month(day)) %>%
     group_by(site, area, mowing, month) %>%
     summarise(frost_days = length(temp_min))
   
   frost_days %>%
     ggplot(aes(x=site, y= frost_days))+
     geom_boxplot()+
     facet_wrap(~month)
   
  ## area 1 C from site 2 is missing from sum_temp
  # calculate mean dff from mown to control in site 2
  gurgl_temp %>%
    filter(site == 2, area == 1, block == "A") %>%
    select(mowing) %>% unique()
  
  sum_temp %>%
    arrange(site, area) %>%
    filter(site == 2 & area != 1) %>%
    group_by(site, area) %>%
    summarize(diffMC_summer = diff(summer_mean)) %>%
    pull(diffMC_summer) %>% mean()  -> mean_diff_CtoM
  
  sum_temp <- sum_temp %>% 
    ungroup() %>%
    add_row(site = "2", area = "1", mowing = "C", snow_days = 142,
            summer_mean = 13.3-mean_diff_CtoM,
            .before = 10)
  
  (sum_temp <- sum_temp %>% arrange(site, area))
  
  # Schrankogel
  # annual
  schran_sumtemp <- schran_temp %>%
    group_by(logger_ID) %>%
    summarise(ann_mean = mean(temp))
  
  # summer
  schran_sumtemp <- schran_temp %>%
    filter(month == 6 | month == 7) %>%
    group_by(logger_ID) %>%
    summarise(summer_mean = mean(temp)) %>%
    full_join(schran_sumtemp)
  
  # number of snow days  
  schran_sumtemp <- schran_temp %>%
    group_by(logger_ID, day) %>%
    summarise(temp_mean = mean(temp),
              temp_min = min(temp),
              temp_max = max(temp)) %>%
    filter(temp_min > -1 &
             temp_max < 1) %>%
    group_by(logger_ID) %>%
    summarise(snow_days = length(temp_mean)) %>%
    full_join(schran_sumtemp)
  
  
  ## day of snowmelt
  schran_sumtemp <- schran_temp %>%
    group_by(logger_ID, day) %>%
    summarise(temp_mean = mean(temp),
              temp_min = min(temp),
              temp_max = max(temp)) %>%
    filter(#temp_mean < 1 & temp_mean > -1 &
      temp_min > -1 &
        temp_max < 1) %>%
    group_by(logger_ID) %>%
    summarise(last_snowday = max(day)) %>%
    right_join(schran_sumtemp)
  
  
  
#### exploratory plots ####
# seeds per species
  all_seeds %>% 
    ggplot(aes(x= reorder(taxon, number_seedlings), y= number_seedlings)) + 
    geom_boxplot() +
    coord_flip()# +
  #  geom_point()

# seeds per site
  tot_seeds %>%
    ggplot(aes(mowing, tot_seeds)) +
    geom_boxplot() +
    facet_grid(~site)
  #pivot_wider(names_from = mowing) %>% 

# mean total seeds per plot
  all_seeds %>% 
    group_by(site, area, mowing) %>%
    summarise(value = sum(number_seedlings)) %>% #group by plots and summarize every plot
    ggplot(aes(mowing, value)) +
    geom_boxplot() +
    facet_grid(~site) 

## mean temperatures
  gurgl_temp %>%
    ggplot(aes(x=mowing, y= temp))+
    geom_boxplot()+
    facet_grid(~site)
  sum_temp %>%
    ggplot(aes(x=mowing, y= summer_mean))+
    geom_boxplot()+
    facet_grid(~site)
  
## temp line through year
  gurgl_temp %>%
    ggplot(aes(x = day, y = temp,
           col = area))+
    geom_point()+
    facet_wrap(~site)
  
# plots
sum_temp %>%
  ggplot(aes(x=mowing, y=summer_mean))+
  geom_boxplot()+
  facet_grid(~site)
sum_temp %>%
  ggplot(aes(x=mowing, y=snow_days))+
  geom_boxplot()+
  facet_grid(~site)

#### climate of previous years to compare with Margreiter 2021

Obergurgl_clim <- as_tibble(read.csv("data/Obergurgel Klimastation Daten 2016-Aug22.csv")%>%
  mutate(date = date(time),
         day = as.factor(day(time)),
         month = as.factor(month(time)),
         year = year(time)) )

Obergurgl_clim %>%
  group_by(year, month) %>%
  summarise(air_temp = mean(t),
            erd_temp = mean(erdmin)) %>%
  ggplot(aes(x=month, y=erd_temp))+
  geom_boxplot()+
  facet_wrap(~year)

Obergurgl_clim %>%
      filter(date >= "2017-05-01" &
               date <= "2017-09-01" |
            date >= "2022-05-01" &
               date <= "2022-09-01") %>%
    mutate(frost_yesno = if_else(erdmin < -0, 1, 0))%>%
  group_by(year) %>%
  summarise(air_temp = mean(t),
            air_temp_min = mean(tmin),
            erd_temp = mean(erdmin),
            frost_days = sum(frost_yesno))
  
  
  
