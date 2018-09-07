library(tidyverse)

# Read in some data. 
network <- read_table2('southernwomen.txt',skip = 2, col_names = FALSE) %>% 
  select(woman = X1, event = X2)

num_women <- unique(network$woman) %>% length()
num_events <- unique(network$event) %>% length()

# convert to all unique IDs
network_edges <- network %>% 
  mutate(
    from = woman,
    to = event + num_women,
  ) %>% 
  select(from, to)


write_csv(network_edges, 'network_edges.csv')
