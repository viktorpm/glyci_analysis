### cheks all the used libraries in the project. Removes duplicates 
library(tidyverse)

map(list.files(pattern = ".R"),
    ~ readLines(con = .x) %>% 
      grep(pattern = "library\\(", x = ., value = T)
) %>%
  unlist() %>%
  substr(start = 1,
         stop = regexpr(text = ., pattern = "\\)") %>% # regexpr: in case of multiple matches it only takes the first one
           as.vector()
  ) %>% `[` (!duplicated(.)) -> liblist

 
### loads all libraries
eval(parse(text=liblist))
