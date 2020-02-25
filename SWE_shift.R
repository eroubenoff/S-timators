# Sweden demo
# Trying to use the shift estimator on Swedish data

SWE <- read_table("SWE_mx.txt", skip = 1)

# Some data reformating
SWE <- SWE %>% mutate(Age = if_else(Age== "110+", "110", Age),
               Age = as.numeric(Age),
               Male = if_else(Male == ".", as.character(NA), Male),
               Male = as.numeric(Male),
               Female = if_else(Female == ".", as.character(NA), Female),
               Female = as.numeric(Female),
               Total = if_else(Total == ".", as.character(NA), Total),
               Total = as.numeric(Total))


m <- SWE %>% filter(Year == 1892) %>%
  select(Age, Male) %>%
  rename(age = Age, v = Male) 

f <- SWE %>% filter(Year == 1892) %>%
  select(Age, Female) %>%
  rename(age = Age, v = Female)

procedure(m, f, 60, 100, 10)  # e1 converges around 11 but e4 to 33
procedure(m, f, 80, 100, 10) #e1 converges to 13 but e2 converes to -8?
procedure(m, f, 80, 90, 100) # e1 converges to 5 but e3 converges to 3.84

# Try again with a rabdin
