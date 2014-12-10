author_details <- function(filename_first_authors, baad) {
  first_authors <- read.csv(filename_first_authors, na.strings=c(NA, ''),
                            stringsAsFactors=FALSE, strip.white=TRUE)
  data_authors <- baad$contact
  data_authors <- data_authors[order(last_name(data_authors$name),
                                     data_authors$name),]
  data_authors$email[data_authors$email == ""] <- NA
  
  cols <- c("name", "email", "address")
  all_authors <- rbind(first_authors[cols], data_authors[,cols])
  all_authors <- all_authors[!duplicated(all_authors$name),]
  
  address <- unique(all_authors$address)
  address_table <- data.frame(code=seq_along(address), address=address)
  
  all_authors$address_code <- match(all_authors$address, address)
  
  data.frame(authors=all_authors,
       address_table=address_table)
}
