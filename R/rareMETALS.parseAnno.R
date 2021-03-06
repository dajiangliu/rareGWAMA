#' A function to pass ANNOFULL from SEQMINER;
#'
#' @param x The input of the ANNOFULL from SEQMINER;
#' @return A string with parsed information;
#' @export
rareMETALS.parseAnno <- function(x){
      if (length(x) == 0) {
              return(NA)
          }
        d <- str_match_all(x, "\\([^\\)]+")[[1]]
        if (length(d) == 0) {
                return(NA)
            }
        ret <- d
        for (idx in 1:length(d)) {
                m <- d[idx]
                codon <- str_split(m, ":")[[1]][1]
                codon <- str_replace(codon, pattern = "^\\(", replacement = "")
                codon <- str_split(codon, "->")[[1]]
                codon <- unlist(sapply(str_split(codon, "/"), `[`, 2))
                
                codonNum <- str_split(m, ":")[[1]]
                codonNum <- str_extract(codonNum, "^Codon[0-9]+")
                codonNum <- codonNum[!is.na(codonNum)]
                codonNum <- str_replace(codonNum, "Codon", "")
                                        #codonNum

                ret[idx] <- str_join(codon[1], codonNum, codon[2])
            }
        ret
  }
