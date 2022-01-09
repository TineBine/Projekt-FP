library(ggdag)
library(ggplot2)
theme_set(theme_dag())

#ce je graf nepregleden, poženemo ukaz še enkrat

#1.graf
dagify(K~6,
       6~5+4,
       5~3+2,
       4~2,
       3~1,
       2~1,
       1~Z) %>% ggdag()

#2.graf
dagify(K~8,
       8~7+6,
       7~5,
       6~4+3,
       5~3+2,
       4~2,
       3~1,
       2~1,
       1~Z) %>% ggdag()

#3.graf
dagify(K~10,
       10~9+8,
       9~7+4,
       8~6+5,
       7~3,
       6~3,
       5~2,
       4~1,
       3~1,
       2~1,
       1~Z) %>% ggdag()

