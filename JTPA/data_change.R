library(dplyr);library(tibble);library(tidyr);library(readr)
# sep = "\t" 明确指定使用制表符作为分隔符
my_data <- read.table("JTPA/jtpa_han.tab", header = TRUE, sep = "\t") 
my_data <- my_data %>% 
  select(-recid,-n_hs2_abv) %>% 
  mutate(
    edu = as.numeric(edu >=12), 
    prevearn = as.vector(scale(log(prevearn+1))),
    n_hs2 = as.numeric(n_hs2 >= 13), # 0    1    2    5    9   13   15   16   17   19   28   37   39   40  346 
    # earnings = as.numeric(earnings>median(earnings))
  )
colnames(my_data) <- c('A1', 'Z1', 'Y2', 'L1', 'A0', 'Z0', 'L0')
write.csv(my_data,'JTPA/jpta_han.csv')
