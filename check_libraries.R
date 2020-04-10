### cheks all the used libraries in the project. Removes duplicates 

map(list.files(pattern = ".R"),
    ~ readLines(con = .x) %>% 
      grep(pattern = "library\\(", x = ., value = T)
) %>%
  unlist() %>%
  substr(start = 1,
         stop = regexpr(text = ., pattern = "\\)") %>% # regexpr: in case of multiple matches it only takes the first one
           as.vector()
  ) %>% `[` (!duplicated(.)) 
