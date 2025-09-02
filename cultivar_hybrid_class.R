### Load libraries
library(ggplot2)
library(dplyr)
setwd("D:/C/Desktop/Olive_GO_paper")
df<-read.table("cult_hybrid_class.txt")
classify_ind <- function(HI, IH){
  # Check pure parents first with low IH
  if(HI >= 0.90 && IH <= 0.20) return("Parent1")
  if(HI <= 0.10 && IH <= 0.20) return("Parent2")
  
  # Early generation hybrids
  if(IH >= 0.75&& HI >= 0.45 && HI <= 0.56) return("F1")
  if(HI >= 0.40 && HI <= 0.60 && IH >= 0.35 && IH < 0.75) return("F2")
  if(HI > 0.60 && IH >= 0.35 && IH <= 0.75) return("BC1_to_P1")
  if(HI < 0.40 && IH >= 0.35 && IH <= 0.75) return("BC1_to_P2")
  if(HI > 0.75 && IH < 0.35 && IH < 0.75) return("BC2_to_P1")
  if(HI < 0.25 && IH < 0.25 ) return("BC2_to_P2")
  

  
  # Ambiguous cases (e.g., low IH with HI near parents or unexpected combos)
  return("Ambiguous")
}
df <- df %>%
  rowwise() %>%
  mutate(Class = classify_ind(HI, IH))

print(df)


write.csv(df, "cultivars_hybrid_classes.csv", row.names = FALSE)
