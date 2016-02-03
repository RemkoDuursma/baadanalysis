
get_wordcount <- function(fn){
  
  # gets first numeric string
  getcount <- function(x)as.numeric(str_extract(x, "[0-9]{1,4}"))
  
  # not robust; will find first matched argument
  count <- function(txt)getcount(r[grep(txt,r)[1]])
  
  if(file.exists(fn)){
    
    r <- readLines(fn)
    
    # should do this more cleverly, but annoyingly the subsections are not added to a total
    # for the section. So more work to make this general
    
    l <- list()
    l$Abstract <- count("Summary")
    l$Introduction <- count("Introduction")
    l$Methods <- count("Data") + count("Data analysis")
    l$Results <- count("Results")
    l$Discussion <- count("Discussion")
    
    l$Total <- count("Words in text")
    
    return(l)
    
  } else {
    return(list(NA))
  }
}
