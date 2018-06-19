#!/usr/bin/env Rscript

# TODO fix final legend order
# TODO plot using hours, but still label days
# TODO put this in the main index.Rmd?
# TODO https://cran.r-project.org/web/packages/cowplot/vignettes/shared_legends.html

suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(tidyr))

#####################
# command-line args #
#####################

read_text <- function(fileName)
  # from https://stackoverflow.com/questions/9068397
  readChar(fileName, file.info(fileName)$size)

args <- docopt(read_text('@usageTxt@'))
od750in <- file.path(getwd(), args[['od750']])
tlistin <- file.path(getwd(), args[['colors']])
plotout <- file.path(getwd(), args[['plots']])

# treatments are guessed from the data, but this sets the factor ordering
# colors is a named vector used to assign colors
treatments <- read.csv(tlistin, header=TRUE, colClasses=c('factor', 'character'), sep='\t')
colors <- treatments$color
names(colors) <- treatments$treatment

find_treatment <- function(treatments, culture) {
  # which treatment does a culture belong to?
  matches <- sapply(treatments, grepl, culture) %>% which %>% names
  # stopifnot(length(matches) == 1)
  return(matches)
}

#########################
# shared graph elements #
#########################

fix_treatment <- function(df) {
  df$treatment <- ordered(df$treatment, levels=treatments$treatment)
  stopifnot(!any(is.na(df$treatment)))
  return(df)
}

# remove_wt <- function(df)
#   subset(df, treatment != 'WT' & treatment != 'KmR') %>% fix_treatment

my_x_scale <- scale_x_continuous(breaks=seq(0,20,by=1))

my_error_bars <- geom_errorbar(size=0.8, width=1.5,
                               position = position_dodge(width=0.5))

# cols <- c('LL'='darkgreen', 'NL'='limegreen', 'HL'='yellowgreen',
#           'WT'='black', 'KmR'='darkgrey')
my_colors  <- scale_colour_manual("Treatment", values=colors)

combine_errors <- function(err1, err2)
  sqrt(err1^2 + err2^2)

#######################
# raw OD measurements #
#######################

od750 <- read.csv(od750in, sep='\t') %>% tbl_df %>% mutate(day=hour/24)
od750 <- od750[!is.na(od750$od750),]
tmp <- sapply(od750$culture, function(c) find_treatment(names(colors), c))
od750$treatment <- tmp

# summarize mean +/- sd per treatment
od750_summary <- od750 %>%
  group_by(treatment, dilution, hour, day) %>%
  summarize(od_mean=mean(od750), od_sd=sd(od750))

plot1 <- od750_summary %>% fix_treatment %>%
  ggplot(aes(x=day, y=od_mean, ymin=od_mean-od_sd, ymax=od_mean+od_sd,
             group=interaction(treatment, dilution), color=treatment)) +
  geom_point(size=2) +
  geom_line(size=0.8) +
  my_x_scale +
  my_error_bars +
  scale_y_continuous(breaks=seq(0, 5, by=0.025)) +
  ylab(bquote('OD'[750])) +
  my_colors

##################################
# doublings within each dilution #
##################################

doublings <- od750 %>%
  group_by(treatment, culture, dilution) %>%
  arrange(hour) %>%
  summarize(initial_od=first(od750), final_od=last(od750), doublings=log2(final_od/initial_od))

doublings_by_dilution <- doublings %>%
  group_by(treatment, dilution) %>%
  summarize(dbl_mean=mean(doublings), dbl_sd=sd(doublings)) 

doublings_total <- doublings_by_dilution %>%
  select(-dbl_sd) %>%
  mutate(dilution=paste0("d", dilution)) %>%
  spread(dilution, dbl_mean) %>%
  mutate(dtot=sum(d1,d2,d3))
rm(doublings, doublings_by_dilution)

doublings_all <- od750 %>%
  group_by(treatment, culture, dilution) %>%
  arrange(hour) %>%
  mutate(od_after_dilution=first(od750),
         doublings_since_dilution=log2(od750/od_after_dilution)) %>% ungroup

print(as.data.frame(doublings_all))

# summarize mean +/- sd per treatment
doublings_all_summary <- doublings_all %>%
  group_by(treatment, dilution, hour, day) %>%
  summarize(dbl_mean=mean(doublings_since_dilution),
            dbl_sd=sd(doublings_since_dilution))

plot2 <- doublings_all_summary %>% fix_treatment %>%
  ggplot(aes(x=day, y=dbl_mean, ymin=dbl_mean-dbl_sd, ymax=dbl_mean+dbl_sd,
             group=interaction(treatment, dilution), color=treatment)) +
    geom_hline(yintercept=0, color='grey', alpha=0.75, linetype='dashed') +
    geom_point(size=2) +
    geom_line(size=0.8) +
    my_error_bars +
    my_x_scale +
    my_colors +
    scale_y_continuous(breaks=seq(-5,10, by=1)) +
    xlab('Total days growth') +
    ylab('Doublings since last dilution')

###################
# total doublings #
###################

# total generations during dilution 1
t1 <- doublings_all %>% filter(dilution == 1) %>%
  mutate(doublings=doublings_since_dilution) %>%
  select(treatment, culture, dilution, hour, day, doublings)

# total generations during dilution 2
t1last <- t1 %>%
  filter(hour == max(hour)) %>% transmute(culture, d1doublings=doublings)
t2 <- doublings_all %>% filter(dilution == 2) %>% left_join(t1last) %>%
  mutate(doublings=d1doublings + doublings_since_dilution) %>%
  select(treatment, culture, dilution, hour, day, doublings)

# total generations during dilution 3
t2last <- t2 %>%
  filter(hour == max(hour)) %>% transmute(culture, d2doublings=doublings)
t3 <- doublings_all %>% filter(dilution == 3) %>% left_join(t2last) %>%
  mutate(doublings=d2doublings + doublings_since_dilution) %>%
  select(treatment, culture, dilution, hour, day, doublings)

# all total generations
doublings_total <- rbind(t1, t2, t3)
rm(t1, t1last, t2, t2last, t3)

# summarize mean +/- sd per treatment
doublings_total_summary <- doublings_total %>%
  group_by(treatment, dilution, hour, day) %>%
  summarize(dbl_mean=mean(doublings), dbl_sd=sd(doublings))

plot3 <- doublings_total_summary %>% fix_treatment %>%
  ggplot(aes(x=day, y=dbl_mean, color=treatment,
             ymin=dbl_mean-dbl_sd, ymax=dbl_mean+dbl_sd)) +
    geom_hline(yintercept=7, color='grey', alpha=0.75, linetype='dashed') +
    geom_point(size=2) +
    geom_line(size=0.8) +
    my_error_bars +
    my_x_scale +
    my_colors +
    scale_y_continuous(breaks=seq(-5,15,by=1)) +
    xlab('Total days growth') + ylab('Doublings total')

###################
# plot everything #
###################

# ggsave(plot1, filename='plot1_meansd.png')
# ggsave(plot2, filename='plot2_meansd.png')
# ggsave(plot3, filename='plot3_meansd.png')

# see https://andyphilips.github.io/blog/2017/04/04/single-legend-for-multiple-plots.html
plot1shared <- plot1 + theme(axis.title.x = element_blank(), legend.position = "none")
plot2shared <- plot2 + theme(axis.title.x = element_blank(), legend.position = "none")
plot3shared <- plot3 + guides(fill="Treatment") + theme(legend.position = c(0.8, 0.3))
plots <- plot_grid(plot1shared, plot2shared, plot3shared, ncol=1, rel_heights = c(1, 1, 2))
ggsave(filename=plotout, plot=plots, width=5, height=10)
