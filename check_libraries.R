### cheks all the used libraries in the project. Removes duplicates 
library(tidyverse)


### list of libraries mentioned in the R scripts of the project
liblist <- map(list.files(pattern = ".Rmd|.R"),
    ~ readLines(con = .x) %>% 
      grep(pattern = "library\\(", x = ., value = T)
) %>%
  unlist() %>%
  substr(start = 9,
         stop = regexpr(text = ., pattern = "\\)") - 1 %>% # regexpr: in case of multiple matches it only takes the first one
           as.vector()
  ) %>% `[` (!duplicated(.))

### list of installed packages
pkgs <- installed.packages() %>% 
  as_tibble() %>% 
  pull(1)

### packages that are needed but not installed
c(liblist[which(!liblist %in% pkgs)])

### installs uninstalled packages
install.packages(c(liblist[which(!liblist %in% pkgs)]))

### loads all libraries

lapply(liblist %>% as.list(), require, character.only = TRUE)

