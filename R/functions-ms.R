author_details <- function(baad, first_authors_file) {
  
  first_authors <- read.csv(first_authors_file, stringsAsFactors=FALSE)
  
  data_authors <- baad$contact
  data_authors <- data_authors[order(last_name(data_authors$name),
                                     data_authors$name),]
  
  cols <- c("name", "address")
  all_authors <- rbind(first_authors[cols], data_authors[,cols])
  all_authors <- all_authors[!duplicated(all_authors$name),]

  all_authors <- subset(all_authors, !name %in% c("Remko A. Duursma","Daniel S. Falster"))
  names(all_authors) <- c("Name","Affiliation")
  
  # split affiliation into multiple rows to avoid nchar > maxchar.
  z <- split(all_authors, 1:nrow(all_authors))
  
  splitaffil <- function(x){
    maxch <- 80
    
    b <- str_trim(strsplit(x[[2]],",")[[1]])
    nc <- cumsum(nchar(b))
    
    txt <- vector("list",10)
    txt[[1]] <- b[1]
    rtxt <- b
    i <- 2
    k <- 1
    for(i in 2:length(rtxt)){
      
      tmp <- paste(paste(txt[[k]], collapse=","), rtxt[i], sep=", ")
      if(nchar(tmp) > maxch){
        k <- k + 1  
      } else {
        txt[[k]] <- c(txt[[k]], rtxt[i])
      }
      
    }
    txt <- txt[!sapply(txt,is.null)]
    if(length(txt) > 1){
      txt <- c(sapply(txt[1:(length(txt)-1)], function(x)paste0(paste(x,collapse=", "),",")),
               paste(txt[[length(txt)]],collapse=", "))
    } else {
      txt <- paste(txt[[1]],collapse=", ")
    }
    
    z <- data.frame(Name=x$Name, Affiliation=txt, stringsAsFactors=FALSE)
    if(nrow(z) > 1)z$Name[2:nrow(z)] <- ""
    return(z)
  }
  l <- lapply(z, splitaffil)
  all_authors <- do.call(rbind,l)
  
  rownames(all_authors) <- as.character(1:nrow(all_authors))
return(all_authors)
}

last_name <- function(author) {
  sapply(author, function(x)
    last(strsplit(x, ' ', useBytes=TRUE)[[1]]))
}

last <- function(x) {
  x[[length(x)]]
}


